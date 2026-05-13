"""
Gaussian vs Linear Smoothing Comparison: v2(pT) for 10x10x10
  - STAR BES-I experimental data
  - 10x10x10 local density (linear block smoothing)
  - 10x10x10 local density (Gaussian kernel smoothing)

Memory optimized: processes ROOT files in chunks to avoid OOM.
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
ROOT_LIN = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/local_density_approach/ana/ampt_localdensity.root'
ROOT_GAU = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/gaussian_correction/ana/ampt_h1.0_200k.root'
STAR_DIR = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/HEPData-ins1395151-v2-csv'
OUTDIR   = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/publication_plots'
os.makedirs(OUTDIR, exist_ok=True)

PID = dict(pip=211, pim=-211, kp=321, km=-321, pr=2212, pbar=-2212)
pt_bins = np.linspace(0.1, 2.4, 24)
pt_centers = 0.5 * (pt_bins[:-1] + pt_bins[1:])

species = [
    ([PID['pip'], PID['pim']], 'Pions',   'Table110.csv', r'STAR $\pi^+$ BES-I'),
    ([PID['kp'],  PID['km']],  'Kaons',   'Table111.csv', r'STAR $K^-$ BES-I'),
    ([PID['pr'],  PID['pbar']],'Protons', 'Table107.csv', r'STAR $p$ BES-I'),
]

# ── Process ROOT file in chunks ──────────────────────────────────────────────
def process_root_in_chunks(filepath, label):
    t0 = time.time()
    
    # Store sums for each species/pt-bin to compute v2 later
    # Format: sums[label] = [sum_v2, sum_v2_sq, count]
    sums = {sp_lbl: [np.zeros(len(pt_centers)), np.zeros(len(pt_centers)), np.zeros(len(pt_centers))] 
            for _, sp_lbl, _, _ in species}
            
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
        
        # Avoid warnings by only computing log where valid
        y[valid] = 0.5 * np.log((e[valid] + pz[valid]) / denom[valid])
        
        pT  = np.sqrt(px**2 + py**2)
        phi = np.arctan2(py, px)
        v2  = np.cos(2 * phi)
        
        for target_pids, sp_lbl, _, _ in species:
            mask = np.isin(pid, target_pids) & (np.abs(y) < 1.0)
            
            pT_sub = pT[mask]
            v2_sub = v2[mask]
            
            # Digitize pT to find which bin each particle falls into (1-indexed)
            indices = np.digitize(pT_sub, pt_bins)
            
            for i in range(1, len(pt_bins)):
                bin_mask = (indices == i)
                bin_v2 = v2_sub[bin_mask]
                
                if len(bin_v2) > 0:
                    sums[sp_lbl][0][i-1] += np.sum(bin_v2)
                    sums[sp_lbl][1][i-1] += np.sum(bin_v2**2)
                    sums[sp_lbl][2][i-1] += len(bin_v2)
                    
    print(f"  {label}: {total_particles} particles processed ({time.time()-t0:.1f}s)")
    
    # Compute final mean and error for each species
    results = {}
    for sp_lbl in sums:
        s_v2, s_v2sq, count = sums[sp_lbl]
        
        mean_v2 = np.full_like(pt_centers, np.nan)
        err_v2 = np.full_like(pt_centers, np.nan)
        
        valid = count > 10
        if np.any(valid):
            mean_v2[valid] = s_v2[valid] / count[valid]
            # Variance = E[X^2] - (E[X])^2
            var_v2 = (s_v2sq[valid] / count[valid]) - (mean_v2[valid])**2
            # Handle float imprecision
            var_v2 = np.maximum(var_v2, 0)
            err_v2[valid] = np.sqrt(var_v2) / np.sqrt(count[valid])
            
        results[sp_lbl] = (mean_v2, err_v2)
        
    return results

# ── Load STAR BES-I data ─────────────────────────────────────────────────────
def load_star_v2(filename):
    filepath = os.path.join(STAR_DIR, filename)
    pt_list, v2_list, err_list = [], [], []
    if not os.path.exists(filepath):
        return None, None, None
        
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

# ── Load data ────────────────────────────────────────────────────────────────
print("Processing ROOT files iteratively to save memory...")
res_lin = process_root_in_chunks(ROOT_LIN, "Linear/Block (200k)")
res_gau = process_root_in_chunks(ROOT_GAU, "Gaussian (200k, h=1.0)")

# ── Plot ─────────────────────────────────────────────────────────────────────
print("Generating plot...")
fig, axs = plt.subplots(1, 3, figsize=(19, 6))

for ax, (_, label, star_file, star_label) in zip(axs, species):

    # STAR BES-I data
    s_pt, s_v2, s_err = load_star_v2(star_file)
    if s_pt is not None:
        ax.errorbar(s_pt, s_v2, yerr=s_err, fmt='o',
                    mfc='none', mec='black', ecolor='black',
                    markersize=8, capsize=3, elinewidth=1.0,
                    label=star_label, zorder=10)
    else:
        print(f"  STAR data missing for {label}")

    # Linear (Block) Smoothing
    v2m_lin, v2e_lin = res_lin[label]
    ax.errorbar(pt_centers, v2m_lin, yerr=v2e_lin,
                color='#2196F3', marker='s', ls='none',
                markersize=6, capsize=2, elinewidth=0.8,
                label=r'Linear Extrapolation (Block)', zorder=8)

    # Gaussian Smoothing
    v2m_gau, v2e_gau = res_gau[label]
    ax.errorbar(pt_centers, v2m_gau, yerr=v2e_gau,
                color='#E91E63', marker='D', ls='none',
                markersize=6, capsize=2, elinewidth=0.8,
                label=r'Gaussian Kernel (h=1.0)', zorder=9)

    # Aesthetics
    ax.axhline(0, color='gray', ls='--', lw=0.8)
    ax.set_xlabel(r'$p_T$ (GeV/c)')
    ax.set_ylabel(r'$v_2 = \langle\cos 2\phi\rangle$')
    ax.set_title(f'$v_2(p_T)$ — {label} (10-40% Centrality)')
    ax.set_xlim(0.0, 2.4)
    ax.set_ylim(-0.02, 0.20)
    
    ax.legend(fontsize=9, loc='upper left')
    ax.tick_params(direction='in', top=True, right=True)

fig.suptitle(
    r'Smoothing Kernel Comparison: $v_2(p_T)$ — $10\times10\times10$ grid vs STAR BES-I @ 7.7 GeV',
    fontsize=15, fontweight='bold', y=1.01
)
fig.tight_layout()
outpath = f'{OUTDIR}/flow_v1_v2_comparison.png'
fig.savefig(outpath, dpi=200, bbox_inches='tight')
print(f"\n✓ Saved: {outpath}")
