"""
Grid Resolution Comparison: v2(pT) for 10x10x10 vs 20x20x20
  - STAR BES-I experimental data
  - 10x10x10 local density (iqmc=2)
  - 20x20x20 local density (iqmc=2)
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import uproot, csv, os, time

# ── Style ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size':   12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'legend.fontsize': 10,
    'lines.linewidth': 1.8,
    'lines.markersize': 6,
    'figure.facecolor': 'white',
})

# ── Paths ────────────────────────────────────────────────────────────────────
ROOT_10 = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/local_density_approach/ana/ampt_localdensity.root'
ROOT_20 = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/kekcc_20x20x20_production/ana/ampt_localdensity.root'
STAR_DIR = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/HEPData-ins1395151-v2-csv'
OUTDIR   = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/publication_plots'
os.makedirs(OUTDIR, exist_ok=True)

PID = dict(pip=211, pim=-211, kp=321, km=-321, pr=2212, pbar=-2212)

# ── Load ROOT file via uproot ────────────────────────────────────────────────
def load_root(filepath, label):
    t0 = time.time()
    with uproot.open(filepath) as f:
        tree = f["ampt"]
        d = tree.arrays(["pid", "px", "py", "pz", "mass"], library="np")
    pid = d["pid"].astype(np.int32)
    px  = d["px"].astype(np.float64)
    py  = d["py"].astype(np.float64)
    pz  = d["pz"].astype(np.float64)
    m   = d["mass"].astype(np.float64)
    e   = np.sqrt(px**2 + py**2 + pz**2 + m**2)
    denom = e - pz
    valid = (denom > 1e-9) & (e > 1e-9)
    y = np.full_like(e, np.nan)
    y[valid] = 0.5 * np.log((e[valid] + pz[valid]) / denom[valid])
    pT  = np.sqrt(px**2 + py**2)
    phi = np.arctan2(py, px)
    v2  = np.cos(2 * phi)
    print(f"  {label}: {len(pid)} particles ({time.time()-t0:.1f}s)")
    return pid, y, pT, v2

# ── Load STAR BES-I data ─────────────────────────────────────────────────────
def load_star_v2(filename):
    filepath = os.path.join(STAR_DIR, filename)
    pt_list, v2_list, err_list = [], [], []
    with open(filepath) as f:
        for row in csv.reader(f):
            if not row or row[0].startswith('#') or row[0].startswith('PT') or row[0].startswith('$'):
                continue
            try:
                pt_val = float(row[0])
                v2_val = row[1]
                if v2_val == '-': continue
                pt_list.append(pt_val)
                v2_list.append(float(v2_val))
                err_list.append(float(row[2]))
            except: continue
    return np.array(pt_list), np.array(v2_list), np.array(err_list)

# ── Compute v2(pT) ──────────────────────────────────────────────────────────
def compute_v2_vs_pt(pid_arr, y_arr, pT_arr, v2_arr, target_pids, pt_bins):
    mask = np.isin(pid_arr, target_pids) & (np.abs(y_arr) < 1.0)
    pT_sub = pT_arr[mask]
    v2_sub = v2_arr[mask]
    ptc = 0.5 * (pt_bins[:-1] + pt_bins[1:])
    v2_mean, v2_err = [], []
    for plo, phi in zip(pt_bins[:-1], pt_bins[1:]):
        in_bin = (pT_sub >= plo) & (pT_sub < phi)
        n = np.sum(in_bin)
        if n > 10:
            v2_mean.append(np.mean(v2_sub[in_bin]))
            v2_err.append(np.std(v2_sub[in_bin]) / np.sqrt(n))
        else:
            v2_mean.append(np.nan)
            v2_err.append(np.nan)
    return ptc, np.array(v2_mean), np.array(v2_err)

# ── Load data ────────────────────────────────────────────────────────────────
print("Loading ROOT files...")
pid10, y10, pT10, v2_10 = load_root(ROOT_10, "10×10×10")
pid20, y20, pT20, v2_20 = load_root(ROOT_20, "20×20×20")

pt_bins = np.linspace(0.1, 2.0, 22)

# ── Define species ───────────────────────────────────────────────────────────
species = [
    ([PID['pip'], PID['pim']], 'Pions',   'Table110.csv', r'STAR $\pi^+$ BES-I'),
    ([PID['kp'],  PID['km']],  'Kaons',   'Table111.csv', r'STAR $K^-$ BES-I'),
    ([PID['pr'],  PID['pbar']],'Protons', 'Table107.csv', r'STAR $p$ BES-I'),
]

# ── Plot ─────────────────────────────────────────────────────────────────────
print("Generating plot...")
fig, axs = plt.subplots(1, 3, figsize=(19, 6))

for ax, (pids_target, label, star_file, star_label) in zip(axs, species):

    # STAR BES-I data
    try:
        s_pt, s_v2, s_err = load_star_v2(star_file)
        ax.errorbar(s_pt, s_v2, yerr=s_err, fmt='o',
                    mfc='none', mec='black', ecolor='black',
                    markersize=8, capsize=3, elinewidth=1.0,
                    label=star_label, zorder=10)
    except Exception as e:
        print(f"  STAR data missing for {label}: {e}")

    # 10x10x10 grid
    ptc, v2m, v2e = compute_v2_vs_pt(pid10, y10, pT10, v2_10, pids_target, pt_bins)
    ax.errorbar(ptc, v2m, yerr=v2e,
                color='#2196F3', marker='s', ls='none',
                markersize=6, capsize=2, elinewidth=0.8,
                label=r'Local Density — $10\times10\times10$ grid', zorder=8)

    # 20x20x20 grid
    ptc, v2m, v2e = compute_v2_vs_pt(pid20, y20, pT20, v2_20, pids_target, pt_bins)
    ax.errorbar(ptc, v2m, yerr=v2e,
                color='#E91E63', marker='D', ls='none',
                markersize=6, capsize=2, elinewidth=0.8,
                label=r'Local Density — $20\times20\times20$ grid', zorder=9)

    # Aesthetics
    ax.axhline(0, color='gray', ls='--', lw=0.8)
    ax.set_xlabel(r'$p_T$ (GeV/c)')
    ax.set_ylabel(r'$v_2 = \langle\cos 2\phi\rangle$')
    ax.set_title(f'$v_2(p_T)$ — {label} (10-40% Centrality)')
    ax.set_xlim(0.1, 2.1)
    ax.legend(fontsize=9, loc='upper left')
    ax.tick_params(direction='in', top=True, right=True)

fig.suptitle(
    r'Grid Resolution Effect on $v_2(p_T)$ — AMPT-CQMF Local Density vs STAR BES-I @ 7.7 GeV',
    fontsize=15, fontweight='bold', y=1.01
)
fig.tight_layout()
outpath = f'{OUTDIR}/grid_resolution_v2_comparison.png'
fig.savefig(outpath, dpi=200, bbox_inches='tight')
print(f"\n✓ Saved: {outpath}")
