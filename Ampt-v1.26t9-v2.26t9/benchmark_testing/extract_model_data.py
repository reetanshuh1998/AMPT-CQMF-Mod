import numpy as np
import uproot
import time
import os
import json

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9'
OUT_FILE = os.path.join(BASE_DIR, 'benchmark_testing', 'model_v2_data.json')

MODELS = {
    'M1_Default': os.path.join(BASE_DIR, 'local_density_approach', 'ana', 'ampt_default.root'),
    'M2_Fixed_rho1': os.path.join(BASE_DIR, 'local_density_approach', 'ana', 'ampt_fixed_rho1.root'),
    'M3_Fixed_rho2': os.path.join(BASE_DIR, 'local_density_approach', 'ana', 'ampt_fixed_rho2.root'),
    'M4_Fixed_rho3': os.path.join(BASE_DIR, 'local_density_approach', 'ana', 'ampt_fixed_rho3.root'),
    'M5_Linear': os.path.join(BASE_DIR, 'local_density_approach', 'ana', 'ampt_localdensity.root'),
    'M6_Gaussian': os.path.join(BASE_DIR, 'gaussian_correction', 'ana', 'ampt_h1.0_200k.root'),
}

PID = dict(pip=211, pim=-211, kp=321, km=-321, pr=2212, pbar=-2212)
pt_bins = np.linspace(0.05, 2.45, 25)
pt_centers = 0.5 * (pt_bins[:-1] + pt_bins[1:])

species = [
    ([PID['pip'], PID['pim']], 'Pions'),
    ([PID['kp'],  PID['km']],  'Kaons'),
    ([PID['pr']],              'Protons'),
    ([PID['pbar']],            'Antiprotons'),
]

# ── Process ROOT file in chunks ──────────────────────────────────────────────
def process_root(filepath, label):
    print(f"Processing {label}...")
    t0 = time.time()
    sums = {sp_lbl: [np.zeros(len(pt_centers)), np.zeros(len(pt_centers)), np.zeros(len(pt_centers))] 
            for _, sp_lbl in species}
            
    total_particles = 0
    try:
        for d in uproot.iterate(f"{filepath}:ampt", ["pid", "px", "py", "pz", "mass"], library="np", step_size="100 MB"):
            pid = d["pid"].astype(np.int32)
            px  = d["px"].astype(np.float64)
            py  = d["py"].astype(np.float64)
            pz  = d["pz"].astype(np.float64)
            m   = d["mass"].astype(np.float64)
            total_particles += len(pid)
            
            e   = np.sqrt(px**2 + py**2 + pz**2 + m**2)
            denom = e - pz
            valid = (denom > 1e-9) & (e > 1e-9)
            y = np.full_like(e, np.nan)
            y[valid] = 0.5 * np.log((e[valid] + pz[valid]) / denom[valid])
            
            pT  = np.sqrt(px**2 + py**2)
            phi = np.arctan2(py, px)
            v2  = np.cos(2 * phi)
            
            for target_pids, sp_lbl in species:
                mask = np.isin(pid, target_pids) & (np.abs(y) < 1.0)
                pT_sub, v2_sub = pT[mask], v2[mask]
                indices = np.digitize(pT_sub, pt_bins)
                for i in range(1, len(pt_bins)):
                    bin_mask = (indices == i)
                    bin_v2 = v2_sub[bin_mask]
                    if len(bin_v2) > 0:
                        sums[sp_lbl][0][i-1] += np.sum(bin_v2)
                        sums[sp_lbl][1][i-1] += np.sum(bin_v2**2)
                        sums[sp_lbl][2][i-1] += len(bin_v2)
    except Exception as e:
        print(f"  Error processing {filepath}: {e}")
        return None
                    
    print(f"  Processed {total_particles} particles in {time.time()-t0:.1f}s")
    
    results = {}
    for sp_lbl in sums:
        s_v2, s_v2sq, count = sums[sp_lbl]
        mean_v2 = np.full_like(pt_centers, np.nan)
        err_v2 = np.full_like(pt_centers, np.nan)
        valid = count > 10
        if np.any(valid):
            mean_v2[valid] = s_v2[valid] / count[valid]
            var_v2 = np.maximum((s_v2sq[valid] / count[valid]) - (mean_v2[valid])**2, 0)
            err_v2[valid] = np.sqrt(var_v2) / np.sqrt(count[valid])
            
        results[sp_lbl] = {
            'pt': pt_centers.tolist(),
            'v2': mean_v2.tolist(),
            'err': err_v2.tolist()
        }
    return results

all_data = {}
for name, filepath in MODELS.items():
    if os.path.exists(filepath):
        res = process_root(filepath, name)
        if res:
            all_data[name] = res
    else:
        print(f"WARNING: File not found: {filepath}")

with open(OUT_FILE, 'w') as f:
    json.dump(all_data, f, indent=2)
    
print(f"\nAll model data saved to {OUT_FILE}")
