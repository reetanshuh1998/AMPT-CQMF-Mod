import numpy as np
import uproot
import time
import os
import json

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9'
OUT_FILE = os.path.join(BASE_DIR, 'benchmark_testing', 'model_v1_data.json')

MODELS = {
    'M1_Default': os.path.join(BASE_DIR, 'local_density_approach', 'ana', 'ampt_default.root'),
    'M2_Fixed_rho1': os.path.join(BASE_DIR, 'local_density_approach', 'ana', 'ampt_fixed_rho1.root'),
    'M3_Fixed_rho2': os.path.join(BASE_DIR, 'local_density_approach', 'ana', 'ampt_fixed_rho2.root'),
    'M4_Fixed_rho3': os.path.join(BASE_DIR, 'local_density_approach', 'ana', 'ampt_fixed_rho3.root'),
    'M5_Linear': os.path.join(BASE_DIR, 'local_density_approach', 'ana', 'ampt_localdensity.root'),
    'M6_Gaussian': os.path.join(BASE_DIR, 'gaussian_correction', 'ana', 'ampt_h1.0_200k.root'),
}

PID = dict(pip=211, pim=-211, kp=321, km=-321, pr=2212, pbar=-2212)
# STAR data is typically plotted from y = -0.8 to +0.8 or so. Let's use bins from -1.0 to 1.0
y_bins = np.linspace(-1.0, 1.0, 21)
y_centers = 0.5 * (y_bins[:-1] + y_bins[1:])

species = [
    ([PID['pip']], 'pip'), ([PID['pim']], 'pim'),
    ([PID['kp']],  'kp'),  ([PID['km']],  'km'),
    ([PID['pr']],  'p'),   ([PID['pbar']],'pbar'),
]

# ── Process ROOT file in chunks ──────────────────────────────────────────────
def process_root(filepath, label):
    print(f"Processing {label}...")
    t0 = time.time()
    sums = {sp_lbl: [np.zeros(len(y_centers)), np.zeros(len(y_centers)), np.zeros(len(y_centers))] 
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
            # STAR directed flow is typically computed w.r.t reaction plane. In AMPT, the reaction plane is at phi=0 (or px axis).
            # So v1 = p_x / p_T or v1 = cos(phi)
            # We use v1 = p_x / p_T for directed flow
            # Wait, usually it's v1 = <p_x / p_T> or <p_x/sqrt(px^2+py^2)> or <cos(phi)>
            # Let's use <cos(phi)> which is standard.
            valid_pt = pT > 1e-9
            v1 = np.full_like(pT, np.nan)
            v1[valid_pt] = px[valid_pt] / pT[valid_pt] # cos(phi)
            
            for target_pids, sp_lbl in species:
                # typically STAR v1 uses a pt cut of 0.4 < pT < 2.0 for protons, and 0.15 < pT < 2.0 for pions
                if 'p' in sp_lbl:
                    mask = np.isin(pid, target_pids) & valid_pt & (pT >= 0.4) & (pT <= 2.0) & ~np.isnan(y)
                else:
                    mask = np.isin(pid, target_pids) & valid_pt & (pT >= 0.15) & (pT <= 2.0) & ~np.isnan(y)
                    
                y_sub, v1_sub = y[mask], v1[mask]
                indices = np.digitize(y_sub, y_bins)
                for i in range(1, len(y_bins)):
                    bin_mask = (indices == i)
                    bin_v1 = v1_sub[bin_mask]
                    if len(bin_v1) > 0:
                        sums[sp_lbl][0][i-1] += np.sum(bin_v1)
                        sums[sp_lbl][1][i-1] += np.sum(bin_v1**2)
                        sums[sp_lbl][2][i-1] += len(bin_v1)
    except Exception as e:
        print(f"  Error processing {filepath}: {e}")
        return None
                    
    print(f"  Processed {total_particles} particles in {time.time()-t0:.1f}s")
    
    results = {}
    for sp_lbl in sums:
        s_v1, s_v1sq, count = sums[sp_lbl]
        mean_v1 = np.full_like(y_centers, np.nan)
        err_v1 = np.full_like(y_centers, np.nan)
        valid = count > 10
        if np.any(valid):
            mean_v1[valid] = s_v1[valid] / count[valid]
            var_v1 = np.maximum((s_v1sq[valid] / count[valid]) - (mean_v1[valid])**2, 0)
            err_v1[valid] = np.sqrt(var_v1) / np.sqrt(count[valid])
            
        results[sp_lbl] = {
            'y': y_centers.tolist(),
            'v1': mean_v1.tolist(),
            'err': err_v1.tolist()
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
    
print(f"\nAll model v1 data saved to {OUT_FILE}")
