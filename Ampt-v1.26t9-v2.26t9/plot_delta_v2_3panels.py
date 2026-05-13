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
    'lines.markersize': 8,
    'figure.facecolor': 'white',
})

# ── Paths ────────────────────────────────────────────────────────────────────
ROOT_LIN = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/local_density_approach/ana/ampt_localdensity.root'
ROOT_GAU = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/gaussian_correction/ana/ampt_h1.0_200k.root'
STAR_DIR = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/HEPData-ins1395151-v2-csv'
STAR_SPLIT_P = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/star_data/v2_splitting_p_pbar_7.7_10_40.csv'

OUTDIR   = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/publication_plots'
os.makedirs(OUTDIR, exist_ok=True)

PID = dict(pip=211, pim=-211, kp=321, km=-321, pr=2212, pbar=-2212)
pt_bins = np.linspace(0.1, 2.0, 20)
pt_centers = 0.5 * (pt_bins[:-1] + pt_bins[1:])

species = [
    ([PID['pip']], 'pip'), ([PID['pim']], 'pim'),
    ([PID['kp']],  'kp'),  ([PID['km']],  'km'),
    ([PID['pr']],  'p'),   ([PID['pbar']],'pbar'),
]

# ── Process ROOT file in chunks ──────────────────────────────────────────────
def process_root_in_chunks(filepath, label):
    sums = {sp_lbl: [np.zeros(len(pt_centers)), np.zeros(len(pt_centers)), np.zeros(len(pt_centers))] 
            for _, sp_lbl in species}
            
    for d in uproot.iterate(f"{filepath}:ampt", ["pid", "px", "py", "pz", "mass"], library="np", step_size="100 MB"):
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
    return results

def load_star_raw(filename):
    filepath = os.path.join(STAR_DIR, filename)
    pt_dict = {}
    with open(filepath) as f:
        for row in csv.reader(f):
            if not row or row[0].startswith('#') or row[0].startswith('PT') or row[0].startswith('$'): continue
            try:
                pt_val = round(float(row[0]), 3)
                v2_val = row[1]
                if v2_val == '-': continue
                pt_dict[pt_val] = (float(v2_val), float(row[2]))
            except: continue
    return pt_dict

def get_star_splitting(f_plus, f_minus):
    dict_p = load_star_raw(f_plus)
    dict_m = load_star_raw(f_minus)
    pts, dv2s, errs = [], [], []
    for pt in sorted(dict_p.keys()):
        if pt in dict_m:
            dv2 = dict_p[pt][0] - dict_m[pt][0]
            err = np.sqrt(dict_p[pt][1]**2 + dict_m[pt][1]**2)
            pts.append(pt)
            dv2s.append(dv2)
            errs.append(err)
    return np.array(pts), np.array(dv2s), np.array(errs)

def load_star_splitting_proton(filepath):
    data = np.genfromtxt(filepath, delimiter=',', skip_header=1)
    return data[:, 0], data[:, 1], data[:, 4]

# ── Load data ────────────────────────────────────────────────────────────────
print("Processing ROOT files...")
res_lin = process_root_in_chunks(ROOT_LIN, "Linear")
res_gau = process_root_in_chunks(ROOT_GAU, "Gaussian")

# STAR Splittings
s_pt_pi, s_dv2_pi, s_err_pi = get_star_splitting('Table110.csv', 'Table109.csv')
s_pt_k,  s_dv2_k,  s_err_k  = get_star_splitting('Table112.csv', 'Table111.csv')
s_pt_p,  s_dv2_p,  s_err_p  = load_star_splitting_proton(STAR_SPLIT_P)

# ── Plot ─────────────────────────────────────────────────────────────────────
print("Generating plot...")
fig, axs = plt.subplots(1, 3, figsize=(18, 5))

configs = [
    (axs[0], 'pip', 'pim', s_pt_pi, s_dv2_pi, s_err_pi, r'$\pi^+ - \pi^-$', [-0.015, 0.015]),
    (axs[1], 'kp',  'km',  s_pt_k,  s_dv2_k,  s_err_k,  r'$K^+ - K^-$',     [-0.04, 0.04]),
    (axs[2], 'p',   'pbar',s_pt_p,  s_dv2_p,  s_err_p,  r'$p - \bar{p}$',   [-0.02, 0.08]),
]

for ax, sp_p, sp_m, spt, sdv2, serr, title, ylim in configs:
    
    # STAR
    if spt is not None and len(spt) > 0:
        ax.errorbar(spt, sdv2, yerr=serr, fmt='k*', mfc='black', mec='black', ecolor='black',
                    label='STAR BES-I', markersize=14, zorder=10)

    # Linear
    v2p_l, e2p_l = res_lin[sp_p]
    v2m_l, e2m_l = res_lin[sp_m]
    dv2_l = v2p_l - v2m_l
    edv2_l = np.sqrt(e2p_l**2 + e2m_l**2)
    valid_l = ~np.isnan(dv2_l) & (pt_centers <= 1.8)
    ax.errorbar(pt_centers[valid_l], dv2_l[valid_l], yerr=edv2_l[valid_l], 
                color='#2196F3', marker='s', ls='none', # REMOVED LINE
                markersize=7, capsize=2, elinewidth=1.2, label='Linear', zorder=8)

    # Gaussian
    v2p_g, e2p_g = res_gau[sp_p]
    v2m_g, e2m_g = res_gau[sp_m]
    dv2_g = v2p_g - v2m_g
    edv2_g = np.sqrt(e2p_g**2 + e2m_g**2)
    valid_g = ~np.isnan(dv2_g) & (pt_centers <= 1.8)
    ax.errorbar(pt_centers[valid_g], dv2_g[valid_g], yerr=edv2_g[valid_g], 
                color='#E91E63', marker='D', ls='none', # REMOVED LINE
                markersize=7, capsize=2, elinewidth=1.2, label='Gaussian', zorder=9)

    ax.axhline(0, color='gray', ls='--', lw=1.0)
    ax.set_xlabel(r'$p_T$ (GeV/$c$)')
    ax.set_ylabel(r'$\Delta v_2$')
    ax.set_title(title + r' Splitting')
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(ylim)
    ax.legend(fontsize=10, loc='upper left')
    ax.tick_params(direction='in', top=True, right=True)

fig.suptitle(r'Particle-Antiparticle $\Delta v_2$ Splitting at $\sqrt{s_{NN}} = 7.7$ GeV (10-40%)', fontsize=16, fontweight='bold', y=1.02)
fig.tight_layout()
outpath = f'{OUTDIR}/delta_v2_all_species.png'
fig.savefig(outpath, dpi=300, bbox_inches='tight')
print(f"\n✓ Saved: {outpath}")
