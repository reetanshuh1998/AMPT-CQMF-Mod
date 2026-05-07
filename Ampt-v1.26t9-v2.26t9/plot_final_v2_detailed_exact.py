import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, time

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 11,
    'axes.labelsize': 13,
    'axes.titlesize': 14,
    'legend.fontsize': 9,
    'figure.dpi': 200,
    'lines.linewidth': 1.8,
})

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'local_density_approach', 'ana')
OUTDIR = 'publication_plots'
os.makedirs(OUTDIR, exist_ok=True)

CONFIGS = {
    'Default (No Medium)': os.path.join(BASE, 'ampt_default.dat'),
    r'Fixed $\rho=1\rho_0$': os.path.join(BASE, 'ampt_fixed_rho1.dat'),
    r'Fixed $\rho=2\rho_0$': os.path.join(BASE, 'ampt_fixed_rho2.dat'),
    r'Fixed $\rho=3\rho_0$': os.path.join(BASE, 'ampt_fixed_rho3.dat'),
    'Local Density Model': os.path.join(BASE, 'ampt_localdensity.dat'),
}

COLORS = {
    'Default (No Medium)': 'royalblue', 
    r'Fixed $\rho=1\rho_0$': 'darkorange',
    r'Fixed $\rho=2\rho_0$': 'forestgreen',
    r'Fixed $\rho=3\rho_0$': 'firebrick',
    'Local Density Model': 'purple'
}
MARKERS = {
    'Default (No Medium)': 'o', 
    r'Fixed $\rho=1\rho_0$': 's',
    r'Fixed $\rho=2\rho_0$': '^',
    r'Fixed $\rho=3\rho_0$': 'D',
    'Local Density Model': 'P'
}

def parse_ampt_fast_all(filepath, b_min=4.5, b_max=9.3):
    if not os.path.exists(filepath): return 0, None, None, None
    pids, pts, v2s, ys = [], [], [], []
    nevents = 0
    n_left = 0
    current_b = 0.0
    with open(filepath, 'r') as f:
        for line in f:
            c = line.split()
            if not c: continue
            if n_left == 0:
                try:
                    n_left = int(c[2])
                    current_b = float(c[3])
                    if b_min <= current_b < b_max: nevents += 1
                except: pass
            else:
                n_left -= 1
                if b_min <= current_b < b_max:
                    try:
                        pid = int(c[0])
                        if abs(pid) in [211, 321, 2212]:
                            px, py, pz, m = map(float, c[1:5])
                            pt = np.sqrt(px**2 + py**2)
                            p2 = pt**2 + pz**2
                            e = np.sqrt(p2 + m**2)
                            if pt > 0 and e > abs(pz):
                                y = 0.5 * np.log((e + pz) / (e - pz))
                                if abs(y) < 1.0:
                                    v2 = (px**2 - py**2) / (pt**2)
                                    pids.append(pid)
                                    pts.append(pt)
                                    v2s.append(v2)
                    except: pass
    return nevents, np.array(pids), np.array(pts), np.array(v2s)

import uproot

def parse_ampt_uproot(filepath, b_min=4.5, b_max=9.3):
    if not os.path.exists(filepath): return 0, None, None, None
    t0 = time.time()
    try:
        with uproot.open(filepath) as f:
            tree = f["ampt"]
            data = tree.arrays(["b", "pid", "px", "py", "pz", "mass"], library="np")
            
        b = data["b"]
        mask_b = (b >= b_min) & (b < b_max)
        
        pid = data["pid"][mask_b]
        px = data["px"][mask_b]
        py = data["py"][mask_b]
        pz = data["pz"][mask_b]
        mass = data["mass"][mask_b]
        
        pt = np.sqrt(px**2 + py**2)
        e = np.sqrt(pt**2 + pz**2 + mass**2)
        
        mask_kin = (pt > 0) & (e > np.abs(pz))
        y = np.full_like(e, np.nan)
        y[mask_kin] = 0.5 * np.log((e[mask_kin] + pz[mask_kin]) / (e[mask_kin] - pz[mask_kin]))
        
        mask_y = np.abs(y) < 1.0
        
        v2 = np.full_like(e, np.nan)
        v2[mask_kin] = (px[mask_kin]**2 - py[mask_kin]**2) / (pt[mask_kin]**2)
        
        final_mask = mask_kin & mask_y
        
        nevents = len(np.unique(b[mask_b])) # Rough estimate
        print(f"Parsed {os.path.basename(filepath)} via UPROOT in {time.time()-t0:.2f}s")
        return nevents, pid[final_mask], pt[final_mask], v2[final_mask]
    except Exception as e:
        print(f"Uproot failed on {filepath}: {e}")
        return 0, None, None, None

def load_star_csv(filename):
    fpath = os.path.join('star_data', filename)
    if not os.path.exists(fpath): return None, None, None
    data = np.genfromtxt(fpath, delimiter=',', skip_header=1)
    return data[:, 0], data[:, 1], data[:, 4]

def binned_mean_exact(pt_array, v2_array, star_pt_centers, width=0.2, min_counts=5):
    v2_means = []
    v2_errs = []
    for c in star_pt_centers:
        mask = (pt_array >= c - width/2) & (pt_array < c + width/2)
        n = np.sum(mask)
        if n >= min_counts:
            v2_means.append(np.mean(v2_array[mask]))
            v2_errs.append(np.std(v2_array[mask]) / np.sqrt(n))
        else:
            v2_means.append(np.nan)
            v2_errs.append(np.nan)
    return np.array(v2_means), np.array(v2_errs)

# Load STAR Data
s_pt_pi, s_v2_pi, s_err_pi = load_star_csv('v2_pip_7.7_0_80.csv')
s_pt_k,  s_v2_k,  s_err_k  = load_star_csv('v2_kp_7.7_0_80.csv')
s_pt_p,  s_v2_p,  s_err_p  = load_star_csv('v2_proton_7.7_10_40.csv')
s_pt_sp, s_v2_sp, s_err_sp = load_star_csv('v2_splitting_p_pbar_7.7_10_40.csv')

cache = {}
for name, fpath in CONFIGS.items():
    # Automatically switch between formats
    if os.path.exists(fpath.replace('.dat', '.root')):
        fpath = fpath.replace('.dat', '.root')
        nev, pids, pts, v2s = parse_ampt_uproot(fpath, 4.5, 9.3)
    else:
        nev, pids, pts, v2s = parse_ampt_fast_all(fpath, 4.5, 9.3)
        
    if nev > 0:
        cache[name] = {'pid': pids, 'pt': pts, 'v2': v2s}

fig, axs = plt.subplots(2, 2, figsize=(14, 10))

for name in CONFIGS.keys():
    if name not in cache: continue
    d = cache[name]
    
    # 1. Pion
    if s_pt_pi is not None:
        sel = np.isin(d['pid'], [211, -211])
        v2_pi, err_pi = binned_mean_exact(d['pt'][sel], d['v2'][sel], s_pt_pi, width=0.2)
        axs[0,0].errorbar(s_pt_pi, v2_pi, yerr=err_pi, fmt=MARKERS[name], color=COLORS[name], label=name, markersize=6)
        
    # 2. Kaon
    if s_pt_k is not None:
        sel = np.isin(d['pid'], [321, -321])
        v2_k, err_k = binned_mean_exact(d['pt'][sel], d['v2'][sel], s_pt_k, width=0.2)
        axs[0,1].errorbar(s_pt_k, v2_k, yerr=err_k, fmt=MARKERS[name], color=COLORS[name], label=name, markersize=6)
        
    # 3. Proton
    if s_pt_p is not None:
        sel = (d['pid'] == 2212)
        v2_p, err_p = binned_mean_exact(d['pt'][sel], d['v2'][sel], s_pt_p, width=0.2)
        axs[1,0].errorbar(s_pt_p, v2_p, yerr=err_p, fmt=MARKERS[name], color=COLORS[name], label=name, markersize=6)
        
    # 4. Proton-Antiproton Splitting
    if s_pt_sp is not None:
        sel_p = (d['pid'] == 2212)
        sel_pb = (d['pid'] == -2212)
        v2_p, err_p = binned_mean_exact(d['pt'][sel_p], d['v2'][sel_p], s_pt_sp, width=0.2, min_counts=5)
        v2_pb, err_pb = binned_mean_exact(d['pt'][sel_pb], d['v2'][sel_pb], s_pt_sp, width=0.2, min_counts=5)
        delta_v2 = v2_p - v2_pb
        delta_err = np.sqrt(err_p**2 + err_pb**2)
        axs[1,1].errorbar(s_pt_sp, delta_v2, yerr=delta_err, fmt=MARKERS[name], color=COLORS[name], label=name, markersize=7)

# Overlay STAR Data
if s_pt_pi is not None: axs[0,0].errorbar(s_pt_pi, s_v2_pi, yerr=s_err_pi, fmt='k*', label='STAR $\pi^+$ (0-80%)', markersize=12, zorder=10)
if s_pt_k is not None: axs[0,1].errorbar(s_pt_k, s_v2_k, yerr=s_err_k, fmt='k*', label='STAR $K^+$ (0-80%)', markersize=12, zorder=10)
if s_pt_p is not None: axs[1,0].errorbar(s_pt_p, s_v2_p, yerr=s_err_p, fmt='k*', label='STAR $p$ (10-40%)', markersize=12, zorder=10)
if s_pt_sp is not None: axs[1,1].errorbar(s_pt_sp, s_v2_sp, yerr=s_err_sp, fmt='k*', label='STAR Splitting (10-40%)', markersize=12, zorder=10)

titles = [r'Pion $v_2$', r'Kaon $v_2$', r'Proton $v_2$', r'Proton-Antiproton Splitting $\Delta v_2$']
for i, ax in enumerate(axs.flatten()):
    ax.axhline(0, color='black', ls='--', lw=1)
    ax.set_xlabel(r'$p_T$ (GeV/$c$)')
    ax.set_title(titles[i])
    ax.legend(loc='upper left')

axs[0,0].set_ylabel(r'$v_2$')
axs[1,0].set_ylabel(r'$v_2$')
axs[1,1].set_ylabel(r'$\Delta v_2 = v_2(p) - v_2(\bar{p})$')
axs[1,1].set_ylim(-0.15, 0.15)
axs[1,1].set_xlim(0.3, 1.8)

fig.suptitle(r'AMPT-CQMF Evaluated at Exact STAR $p_T$ Bins (7.7 GeV)', fontsize=16, fontweight='bold', y=0.98)
fig.tight_layout()

out_path = os.path.join(OUTDIR, 'exact_pt_detailed_v2.png')
fig.savefig(out_path, bbox_inches='tight')
print(f"\nSUCCESS: Exact pT detailed plot saved to {out_path}")
