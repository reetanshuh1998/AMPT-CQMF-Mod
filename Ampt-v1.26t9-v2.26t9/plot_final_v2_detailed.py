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
    'Local Density Model': os.path.join(BASE, 'ampt_localdensity.dat'),
}

COLORS = {'Default (No Medium)': 'royalblue', 'Local Density Model': 'firebrick'}
MARKERS = {'Default (No Medium)': 'o', 'Local Density Model': 's'}

def parse_ampt_fast_all(filepath, b_min=4.5, b_max=9.3):
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
                        # Keep Pions, Kaons, Protons
                        if abs(pid) in [211, 321, 2212]:
                            pids.append(pid)
                            pxs.append(float(cols[1]))
                            pys.append(float(cols[2]))
                            pzs.append(float(cols[3]))
                            masses.append(float(cols[4]))
                    except: pass
                    
    print(f"Parsed {os.path.basename(filepath)}: {nevents} events in 10-40% centrality ({time.time()-t0:.1f}s)")
    return nevents, np.array(pids), np.array(pxs), np.array(pys), np.array(pzs), np.array(masses)

def binned_mean(x, y, bins, min_counts=10):
    centers = 0.5 * (bins[:-1] + bins[1:])
    means = np.full(len(centers), np.nan)
    errs  = np.full(len(centers), np.nan)
    for i in range(len(bins)-1):
        mask = (x >= bins[i]) & (x < bins[i+1]) & np.isfinite(y)
        n = np.sum(mask)
        if n >= min_counts:
            means[i] = np.mean(y[mask])
            errs[i]  = np.std(y[mask]) / np.sqrt(n)
    return centers, means, errs

def load_star_csv(filename):
    fpath = os.path.join('star_data', filename)
    if not os.path.exists(fpath): return None, None, None
    try:
        data = np.genfromtxt(fpath, delimiter=',', skip_header=1)
        return data[:, 0], data[:, 1], data[:, 4] # pt, v2, stat error
    except: return None, None, None

fig, axs = plt.subplots(2, 2, figsize=(14, 10))
pt_bins = np.linspace(0.2, 2.2, 10)

cache = {}
import uproot
def parse_ampt_uproot(filepath, b_min=4.5, b_max=9.3):
    if not os.path.exists(filepath): return 0, None, None, None, None, None
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
        nevents = len(np.unique(b[mask_b]))
        print(f"Parsed {os.path.basename(filepath)} via UPROOT in {time.time()-t0:.2f}s")
        return nevents, pid, px, py, pz, mass
    except Exception as e:
        print(f"Uproot failed: {e}")
        return 0, None, None, None, None, None

for name, fpath in CONFIGS.items():
    if os.path.exists(fpath.replace('.dat', '.root')):
        fpath = fpath.replace('.dat', '.root')
        nev, pids, pxs, pys, pzs, masses = parse_ampt_uproot(fpath, 4.5, 9.3)
    else:
        nev, pids, pxs, pys, pzs, masses = parse_ampt_fast_all(fpath, 4.5, 9.3)
    if nev > 0:
        pt = np.sqrt(pxs**2 + pys**2)
        e = np.sqrt(pxs**2 + pys**2 + pzs**2 + masses**2)
        with np.errstate(divide='ignore', invalid='ignore'):
            y = 0.5 * np.log((e + pzs) / (e - pzs))
            v2 = (pxs**2 - pys**2) / (pxs**2 + pys**2)
        mask = (np.abs(y) < 1.0) & (pt > 0) & np.isfinite(v2)
        cache[name] = {'pt': pt[mask], 'v2': v2[mask], 'pid': pids[mask]}

# --- Plotting ---
for ci, name in enumerate(CONFIGS.keys()):
    if name not in cache: continue
    d = cache[name]
    
    # Panel 1: Pion v2 (using 0-80% data for comparison though our cut is 10-40%)
    sel_pi = np.isin(d['pid'], [211, -211])
    cen, v2_pi, err_pi = binned_mean(d['pt'][sel_pi], d['v2'][sel_pi], pt_bins, min_counts=50)
    axs[0,0].errorbar(cen, v2_pi, yerr=err_pi, fmt=MARKERS[name], color=COLORS[name], label=name, markersize=6)
    
    # Panel 2: Kaon v2
    sel_k = np.isin(d['pid'], [321, -321])
    cen, v2_k, err_k = binned_mean(d['pt'][sel_k], d['v2'][sel_k], pt_bins, min_counts=50)
    axs[0,1].errorbar(cen, v2_k, yerr=err_k, fmt=MARKERS[name], color=COLORS[name], label=name, markersize=6)
    
    # Panel 3: Proton v2 Absolute
    sel_p = (d['pid'] == 2212)
    cen, v2_p, err_p = binned_mean(d['pt'][sel_p], d['v2'][sel_p], pt_bins, min_counts=20)
    axs[1,0].errorbar(cen, v2_p, yerr=err_p, fmt=MARKERS[name], color=COLORS[name], label=name, markersize=6)
    
    # Panel 4: Proton-Antiproton Splitting
    sel_pb = (d['pid'] == -2212)
    _, v2_pb, err_pb = binned_mean(d['pt'][sel_pb], d['v2'][sel_pb], pt_bins, min_counts=10) # Antiprotons are rare
    
    delta_v2 = v2_p - v2_pb
    delta_err = np.sqrt(err_p**2 + err_pb**2)
    axs[1,1].errorbar(cen, delta_v2, yerr=delta_err, fmt=MARKERS[name], color=COLORS[name], label=name, markersize=7)

# Overlay STAR Data
s_pt, s_v2, s_err = load_star_csv('v2_pip_7.7_0_80.csv')
if s_pt is not None: axs[0,0].errorbar(s_pt, s_v2, yerr=s_err, fmt='k*', label='STAR $\pi^+$ (0-80%)', markersize=10, zorder=10)

s_pt, s_v2, s_err = load_star_csv('v2_kp_7.7_0_80.csv')
if s_pt is not None: axs[0,1].errorbar(s_pt, s_v2, yerr=s_err, fmt='k*', label='STAR $K^+$ (0-80%)', markersize=10, zorder=10)

s_pt, s_v2, s_err = load_star_csv('v2_proton_7.7_10_40.csv')
if s_pt is not None: axs[1,0].errorbar(s_pt, s_v2, yerr=s_err, fmt='k*', label='STAR $p$ (10-40%)', markersize=10, zorder=10)

s_pt, s_dv2, s_err = load_star_csv('v2_splitting_p_pbar_7.7_10_40.csv')
if s_pt is not None: axs[1,1].errorbar(s_pt, s_dv2, yerr=s_err, fmt='k*', label='STAR Splitting (10-40%)', markersize=10, zorder=10)

# Formatting
titles = [r'Pion $v_2(p_T)$', r'Kaon $v_2(p_T)$', r'Proton Absolute $v_2(p_T)$', r'Proton-Antiproton Splitting $\Delta v_2$']
for i, ax in enumerate(axs.flatten()):
    ax.axhline(0, color='black', ls='--', lw=1)
    ax.set_xlabel(r'$p_T$ (GeV/$c$)')
    ax.set_title(titles[i])
    ax.legend(loc='upper left')

axs[0,0].set_ylabel(r'$v_2$')
axs[1,0].set_ylabel(r'$v_2$')
axs[1,1].set_ylabel(r'$\Delta v_2 = v_2(p) - v_2(\bar{p})$')
axs[1,1].set_ylim(-0.02, 0.08)

fig.suptitle(r'AMPT-CQMF vs STAR at 7.7 GeV Au+Au', fontsize=16, fontweight='bold', y=0.98)
fig.tight_layout()

out_path = os.path.join(OUTDIR, 'detailed_validation_v2.png')
fig.savefig(out_path, bbox_inches='tight')
print(f"\nSUCCESS: Detailed validation plot saved to {out_path}")
