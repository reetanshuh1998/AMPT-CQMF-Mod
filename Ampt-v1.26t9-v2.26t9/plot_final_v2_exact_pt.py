import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

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

def load_star_csv(filename):
    fpath = os.path.join('star_data', filename)
    if not os.path.exists(fpath): return None, None, None
    data = np.genfromtxt(fpath, delimiter=',', skip_header=1)
    return data[:, 0], data[:, 1], data[:, 4]

# We will match the AMPT data EXACTLY to the STAR pt points
# STAR Delta v2 points:
s_pt, s_dv2, s_err = load_star_csv('v2_splitting_p_pbar_7.7_10_40.csv')

# Define bin edges so the bin centers exactly match STAR data points
pt_bin_width = 0.3
pt_bins = []
for p in s_pt:
    pt_bins.append([p - pt_bin_width/2, p + pt_bin_width/2])

def get_v2_in_bins(fpath):
    # Parse exactly as before
    pids, pts, v2s, ys = [], [], [], []
    with open(fpath) as f:
        n_left = 0
        current_b = 0.0
        for line in f:
            c = line.split()
            if not c: continue
            if n_left == 0:
                try:
                    n_left = int(c[2])
                    current_b = float(c[3])
                except: pass
            else:
                n_left -= 1
                if 4.5 <= current_b < 9.3:
                    try:
                        pid = int(c[0])
                        if abs(pid) == 2212:
                            px, py, pz, m = map(float, c[1:5])
                            pt = np.sqrt(px**2 + py**2)
                            p2 = px**2 + py**2 + pz**2
                            e = np.sqrt(p2 + m**2)
                            if pt > 0 and e > abs(pz):
                                y = 0.5 * np.log((e + pz) / (e - pz))
                                if abs(y) < 1.0:
                                    v2 = (px**2 - py**2) / (pt**2)
                                    pids.append(pid)
                                    pts.append(pt)
                                    v2s.append(v2)
                    except: pass
    pids, pts, v2s = np.array(pids), np.array(pts), np.array(v2s)
    
    v2_p = []
    v2_pb = []
    err_p = []
    err_pb = []
    
    for (low, high) in pt_bins:
        mask = (pts >= low) & (pts < high)
        m_p = mask & (pids == 2212)
        m_pb = mask & (pids == -2212)
        
        if np.sum(m_p) > 5:
            v2_p.append(np.mean(v2s[m_p]))
            err_p.append(np.std(v2s[m_p]) / np.sqrt(np.sum(m_p)))
        else:
            v2_p.append(np.nan)
            err_p.append(np.nan)
            
        if np.sum(m_pb) > 5:
            v2_pb.append(np.mean(v2s[m_pb]))
            err_pb.append(np.std(v2s[m_pb]) / np.sqrt(np.sum(m_pb)))
        else:
            v2_pb.append(np.nan)
            err_pb.append(np.nan)
            
    return np.array(v2_p), np.array(err_p), np.array(v2_pb), np.array(err_pb)

fig, ax = plt.subplots(figsize=(8, 6))

for name, fpath in CONFIGS.items():
    v2_p, err_p, v2_pb, err_pb = get_v2_in_bins(fpath)
    
    delta_v2 = v2_p - v2_pb
    delta_err = np.sqrt(err_p**2 + err_pb**2)
    
    # Plot with a slight horizontal offset so points don't overlap perfectly
    offset = -0.02 if name == 'Default (No Medium)' else 0.02
    ax.errorbar(s_pt, delta_v2, yerr=delta_err, fmt=MARKERS[name], 
                color=COLORS[name], label=name, markersize=8, capsize=3)

# STAR Delta v2 Data
ax.errorbar(s_pt, s_dv2, yerr=s_err, fmt='k*', label='STAR Au+Au 7.7 GeV (10-40%)', markersize=14, capsize=4, zorder=10)

ax.axhline(0, color='black', ls='--', lw=1)
ax.set_xlabel(r'$p_T$ (GeV/$c$)')
ax.set_ylabel(r'$\Delta v_2 = v_2(p) - v_2(\bar{p})$')
ax.set_title(r'Exact $p_T$ Match: Proton-Antiproton $v_2$ Splitting')
ax.legend(loc='best')

# Zoom into the relevant physics region
ax.set_ylim(-0.15, 0.15)
ax.set_xlim(0.3, 1.8)

fig.tight_layout()
out_path = os.path.join(OUTDIR, 'exact_pt_v2_splitting.png')
fig.savefig(out_path, bbox_inches='tight')
print(f"\nSUCCESS: Exact pT validation plot saved to {out_path}")
