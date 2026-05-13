import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import uproot, csv, os, time

# ── Style ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size':   14,
    'axes.labelsize': 16,
    'axes.titlesize': 16,
    'legend.fontsize': 12,
    'lines.linewidth': 2.0,
    'lines.markersize': 8,
    'figure.facecolor': 'white',
})

# ── Paths ────────────────────────────────────────────────────────────────────
ROOT_LIN = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/local_density_approach/ana/ampt_localdensity.root'
ROOT_GAU = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/gaussian_correction/ana/ampt_h1.0_200k.root'
STAR_SPLIT = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/star_data/v2_splitting_p_pbar_7.7_10_40.csv'

OUTDIR   = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/publication_plots'
os.makedirs(OUTDIR, exist_ok=True)

PID = dict(pr=2212, pbar=-2212)
# We only plot up to ~1.8 to match STAR data range
pt_bins = np.linspace(0.1, 1.9, 10)
pt_centers = 0.5 * (pt_bins[:-1] + pt_bins[1:])

species = [
    ([PID['pr']],   'Protons'),
    ([PID['pbar']], 'Antiprotons'),
]

# ── Process ROOT file in chunks ──────────────────────────────────────────────
def process_root_in_chunks(filepath, label):
    t0 = time.time()
    sums = {sp_lbl: [np.zeros(len(pt_centers)), np.zeros(len(pt_centers)), np.zeros(len(pt_centers))] 
            for _, sp_lbl in species}
            
    total_particles = 0
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
                    
    print(f"  {label}: {total_particles} particles processed ({time.time()-t0:.1f}s)")
    
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
        results[sp_lbl] = (mean_v2, err_v2)
        
        # Print antiproton statistics to help explain the error bars
        if sp_lbl == 'Antiprotons':
            print(f"    Antiproton stats per bin: {count}")
            
    return results

def load_star_splitting(filepath):
    if not os.path.exists(filepath): return None, None, None
    data = np.genfromtxt(filepath, delimiter=',', skip_header=1)
    return data[:, 0], data[:, 1], data[:, 4]

# ── Load data ────────────────────────────────────────────────────────────────
print("Processing ROOT files...")
res_lin = process_root_in_chunks(ROOT_LIN, "Linear/Block (200k)")
res_gau = process_root_in_chunks(ROOT_GAU, "Gaussian (200k, h=1.0)")

# ── Plot ─────────────────────────────────────────────────────────────────────
print("Generating plot...")
fig, ax = plt.subplots(figsize=(8, 6))

# Overlay STAR Data
s_pt_sp, s_v2_sp, s_err_sp = load_star_splitting(STAR_SPLIT)
if s_pt_sp is not None:
    # Plot STAR data with sys error bars if needed, or just stat. We use stat here.
    ax.errorbar(s_pt_sp, s_v2_sp, yerr=s_err_sp, fmt='k*', mfc='black', mec='black', ecolor='black',
                label='STAR BES-I (10-40%)', markersize=14, zorder=10)

# Linear Splitting
v2p_lin, e2p_lin = res_lin['Protons']
v2pb_lin, e2pb_lin = res_lin['Antiprotons']
dv2_lin = v2p_lin - v2pb_lin
edv2_lin = np.sqrt(e2p_lin**2 + e2pb_lin**2)

# Mask out empty bins
valid_lin = ~np.isnan(dv2_lin) & (pt_centers <= 1.8)
ax.errorbar(pt_centers[valid_lin], dv2_lin[valid_lin], yerr=edv2_lin[valid_lin], 
            color='#2196F3', marker='s', ls='-', lw=1.5,
            markersize=8, capsize=3, elinewidth=1.5, label=r'Linear Extrapolation (Block)', zorder=8)

# Gaussian Splitting
v2p_gau, e2p_gau = res_gau['Protons']
v2pb_gau, e2pb_gau = res_gau['Antiprotons']
dv2_gau = v2p_gau - v2pb_gau
edv2_gau = np.sqrt(e2p_gau**2 + e2pb_gau**2)

# Mask out empty bins
valid_gau = ~np.isnan(dv2_gau) & (pt_centers <= 1.8)
ax.errorbar(pt_centers[valid_gau], dv2_gau[valid_gau], yerr=edv2_gau[valid_gau], 
            color='#E91E63', marker='D', ls='-', lw=1.5,
            markersize=8, capsize=3, elinewidth=1.5, label=r'Gaussian Kernel (h=1.0 fm)', zorder=9)

ax.axhline(0, color='gray', ls='--', lw=1.0)
ax.set_xlabel(r'$p_T$ (GeV/$c$)')
ax.set_ylabel(r'$\Delta v_2 = v_2(p) - v_2(\bar{p})$')
ax.set_title(r'Proton-Antiproton Splitting $\Delta v_2$ at $\sqrt{s_{NN}} = 7.7$ GeV')
ax.set_xlim(0.1, 1.9)
ax.set_ylim(-0.02, 0.10)
ax.legend(fontsize=11, loc='upper right')
ax.tick_params(direction='in', top=True, right=True)

fig.tight_layout()
outpath = f'{OUTDIR}/delta_v2_only_comparison.png'
fig.savefig(outpath, dpi=300, bbox_inches='tight')
print(f"\n✓ Saved: {outpath}")
