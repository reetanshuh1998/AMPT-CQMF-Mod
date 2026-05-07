import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, time

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 15,
    'legend.fontsize': 11,
    'figure.dpi': 200,
    'lines.linewidth': 2.0,
})

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'local_density_approach', 'ana')
OUTDIR = 'publication_plots'
os.makedirs(OUTDIR, exist_ok=True)

CONFIGS = {
    'Default (No Medium)': os.path.join(BASE, 'ampt_default.dat'),
    'Local Density Model': os.path.join(BASE, 'ampt_localdensity.dat'),
}

COLORS = {'Default (No Medium)': 'royalblue', 'Local Density Model': 'firebrick'}
MARKERS = {'Default (No Medium)': 'o', 'Local Density Model': 's'}

def parse_ampt_fast_protons(filepath, b_min=4.5, b_max=9.3):
    if not os.path.exists(filepath):
        print(f"Missing {filepath}")
        return 0, None, None, None, None, None
        
    t0 = time.time()
    pids, pxs, pys, pzs, masses = [], [], [], [], []
    nevents = 0
    n_left = 0
    current_b = 0.0
    
    with open(filepath, 'r') as f:
        for line in f:
            cols = line.split()
            if not cols: continue
            if n_left == 0:
                try:
                    n_left = int(cols[2])
                    current_b = float(cols[3])
                    if b_min <= current_b < b_max:
                        nevents += 1
                except: pass
            else:
                n_left -= 1
                if b_min <= current_b < b_max:
                    try:
                        pid = int(cols[0])
                        if abs(pid) == 2212: # Only grab protons
                            pids.append(pid)
                            pxs.append(float(cols[1]))
                            pys.append(float(cols[2]))
                            pzs.append(float(cols[3]))
                            masses.append(float(cols[4]))
                    except: pass
                    
    print(f"Parsed {os.path.basename(filepath)}: {nevents} events in 10-40% centrality ({time.time()-t0:.1f}s)")
    return nevents, np.array(pids), np.array(pxs), np.array(pys), np.array(pzs), np.array(masses)

def binned_mean(x, y, bins):
    centers = 0.5 * (bins[:-1] + bins[1:])
    means = np.full(len(centers), np.nan)
    errs  = np.full(len(centers), np.nan)
    for i in range(len(bins)-1):
        mask = (x >= bins[i]) & (x < bins[i+1]) & np.isfinite(y)
        n = np.sum(mask)
        if n > 50:
            means[i] = np.mean(y[mask])
            errs[i]  = np.std(y[mask]) / np.sqrt(n)
    return centers, means, errs

def load_star_csv(filename):
    fpath = os.path.join('star_data', filename)
    if not os.path.exists(fpath): return None, None, None
    try:
        data = np.genfromtxt(fpath, delimiter=',', skip_header=1)
        return data[:, 0], data[:, 1], data[:, 4]
    except: return None, None, None

fig, axs = plt.subplots(1, 2, figsize=(15, 6))
pt_bins = np.linspace(0.4, 2.0, 8)

import uproot
def parse_ampt_uproot_protons(filepath, b_min=4.5, b_max=9.3):
    if not os.path.exists(filepath): return 0, None, None, None, None, None
    t0 = time.time()
    try:
        with uproot.open(filepath) as f:
            tree = f["ampt"]
            data = tree.arrays(["b", "pid", "px", "py", "pz", "mass"], library="np")
        b = data["b"]
        mask_b = (b >= b_min) & (b < b_max)
        pid = data["pid"][mask_b]
        mask_p = np.abs(pid) == 2212
        nevents = len(np.unique(b[mask_b]))
        print(f"Parsed {os.path.basename(filepath)} via UPROOT in {time.time()-t0:.2f}s")
        return nevents, pid[mask_p], data["px"][mask_b][mask_p], data["py"][mask_b][mask_p], data["pz"][mask_b][mask_p], data["mass"][mask_b][mask_p]
    except Exception as e:
        print(f"Uproot failed: {e}")
        return 0, None, None, None, None, None

for name, fpath in CONFIGS.items():
    if os.path.exists(fpath.replace('.dat', '.root')):
        fpath = fpath.replace('.dat', '.root')
        nev, pids, pxs, pys, pzs, masses = parse_ampt_uproot_protons(fpath, 4.5, 9.3)
    else:
        nev, pids, pxs, pys, pzs, masses = parse_ampt_fast_protons(fpath, 4.5, 9.3)
    if nev == 0: continue
    
    pt = np.sqrt(pxs**2 + pys**2)
    e = np.sqrt(pxs**2 + pys**2 + pzs**2 + masses**2)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        y = 0.5 * np.log((e + pzs) / (e - pzs))
        v2 = (pxs**2 - pys**2) / (pxs**2 + pys**2)
        
    mask = (np.abs(y) < 1.0) & (pt > 0) & np.isfinite(v2)
    
    sel_p = mask & (pids == 2212)
    sel_pb = mask & (pids == -2212)
    
    cen, v2_p, err_p = binned_mean(pt[sel_p], v2[sel_p], pt_bins)
    _, v2_pb, err_pb = binned_mean(pt[sel_pb], v2[sel_pb], pt_bins)
    
    # Right panel: v2
    axs[1].errorbar(cen, v2_p, yerr=err_p, fmt=MARKERS[name]+'-', color=COLORS[name], label=f'{name} (p)', capsize=3)
    axs[1].errorbar(cen, v2_pb, yerr=err_pb, fmt=MARKERS[name]+'--', color=COLORS[name], markerfacecolor='white', label=f'{name} ($\\bar{{p}}$)', capsize=3)
    
    # Left panel: Delta v2
    delta_v2 = v2_p - v2_pb
    delta_err = np.sqrt(err_p**2 + err_pb**2)
    axs[0].errorbar(cen, delta_v2, yerr=delta_err, fmt=MARKERS[name]+'-', color=COLORS[name], label=name, markersize=8, linewidth=2)

# STAR Delta v2 Data
s_pt, s_dv2, s_err = load_star_csv('v2_splitting_p_pbar_7.7_10_40.csv')
if s_pt is not None:
    axs[0].errorbar(s_pt, s_dv2, yerr=s_err, fmt='k*', label='STAR Au+Au 7.7 GeV (10-40%)', markersize=14, capsize=4, zorder=10)

axs[0].axhline(0, color='black', ls='--', lw=1)
axs[0].set_xlabel(r'$p_T$ (GeV/$c$)')
axs[0].set_ylabel(r'$\Delta v_2 = v_2(p) - v_2(\bar{p})$')
axs[0].set_title(r'Proton-Antiproton Elliptic Flow Splitting')
axs[0].legend(loc='upper left')
axs[0].set_ylim(-0.02, 0.08)

axs[1].axhline(0, color='black', ls='--', lw=1)
axs[1].set_xlabel(r'$p_T$ (GeV/$c$)')
axs[1].set_ylabel(r'$v_2$')
axs[1].set_title(r'Absolute $v_2(p_T)$')
axs[1].legend(loc='upper left')

fig.suptitle(r'Validation: Local Density Mean-Field Effect on $v_2$ at 7.7 GeV', fontsize=18, fontweight='bold', y=0.98)
fig.tight_layout()

out_path = os.path.join(OUTDIR, 'final_validation_v2_splitting.png')
fig.savefig(out_path, bbox_inches='tight')
print(f"\nSUCCESS: Validation plot saved to {out_path}")
