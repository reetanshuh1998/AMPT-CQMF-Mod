import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import uproot
import os
import csv
from scipy.stats import linregress

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
STAR_DIR = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/HEPData-ins1277069-v1-csv'
OUTDIR   = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/publication_plots'
os.makedirs(OUTDIR, exist_ok=True)

PID = dict(pip=211, pim=-211, pr=2212, pbar=-2212)
species = [
    ([PID['pip']], 'pip'), ([PID['pim']], 'pim'),
    ([PID['pr']],  'p'),   ([PID['pbar']],'pbar'),
]

# ── Load STAR Data ───────────────────────────────────────────────────────────
def load_star_v1(table_name, target_energy=7.7, particle_str='proton'):
    filepath = os.path.join(STAR_DIR, table_name)
    data = []
    in_target_block = False
    correct_particle = False
    
    if not os.path.exists(filepath): return np.array([]), np.array([]), np.array([])
        
    with open(filepath) as f:
        reader = csv.reader(f)
        for row in reader:
            if not row: continue
            if row[0].startswith('#: SQRT(sNN)'):
                energy = float(row[1])
                in_target_block = np.isclose(energy, target_energy)
                continue
            if row[0].startswith('rapidity'):
                if particle_str == 'proton':
                    correct_particle = ('proton' in row[1] and 'anti' not in row[1])
                else:
                    correct_particle = (particle_str in row[1])
                continue
            if row[0].startswith('#') or row[0].startswith('$'): continue
            if in_target_block and correct_particle:
                try: data.append([float(row[0]), float(row[1]), float(row[2])])
                except ValueError: continue
                    
    data = np.array(data)
    if len(data) == 0: return np.array([]), np.array([]), np.array([])
    return data[:, 0], data[:, 1], data[:, 2]

s_y_p, s_v1_p, s_err_p = load_star_v1('Table3.csv', 7.7, 'proton')
s_y_pbar, s_v1_pbar, s_err_pbar = load_star_v1('Table3.csv', 7.7, 'anti-proton')
s_y_pip, s_v1_pip, s_err_pip = load_star_v1('Table4.csv', 7.7, '\\pi^{+}')
s_y_pim, s_v1_pim, s_err_pim = load_star_v1('Table4.csv', 7.7, '\\pi^{-}')

# Calculate Delta v1
if len(s_y_p) > 0 and len(s_y_pbar) > 0 and np.array_equal(s_y_p, s_y_pbar):
    s_dv1_p = s_v1_p - s_v1_pbar
    s_edv1_p = np.sqrt(s_err_p**2 + s_err_pbar**2)
else:
    s_dv1_p, s_edv1_p = np.array([]), np.array([])

if len(s_y_pip) > 0 and len(s_y_pim) > 0 and np.array_equal(s_y_pip, s_y_pim):
    s_dv1_pi = s_v1_pip - s_v1_pim
    s_edv1_pi = np.sqrt(s_err_pip**2 + s_err_pim**2)
else:
    s_dv1_pi, s_edv1_pi = np.array([]), np.array([])

y_centers_dict = {
    'pip': s_y_pip if len(s_y_pip)>0 else np.array([-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9]),
    'pim': s_y_pim if len(s_y_pim)>0 else np.array([-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9]),
    'p': s_y_p if len(s_y_p)>0 else np.array([-0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7]),
    'pbar': s_y_pbar if len(s_y_pbar)>0 else np.array([-0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7]),
}

def binned_mean_exact(filepath, label):
    sums = {sp_lbl: [np.zeros(len(y_centers_dict[sp_lbl])), np.zeros(len(y_centers_dict[sp_lbl])), np.zeros(len(y_centers_dict[sp_lbl]))] 
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
        valid_pt = pT > 1e-9
        v1 = np.full_like(pT, np.nan)
        v1[valid_pt] = px[valid_pt] / pT[valid_pt]
        
        for target_pids, sp_lbl in species:
            if 'p' in sp_lbl and sp_lbl != 'pip' and sp_lbl != 'pim':
                mask = np.isin(pid, target_pids) & valid_pt & (pT >= 0.4) & (pT <= 2.0) & ~np.isnan(y)
            else:
                mask = np.isin(pid, target_pids) & valid_pt & (pT >= 0.15) & (pT <= 2.0) & ~np.isnan(y)
                
            y_sub, v1_sub = y[mask], v1[mask]
            y_centers = y_centers_dict[sp_lbl]
            
            for i, y_cen in enumerate(y_centers):
                bin_mask = (y_sub >= y_cen - 0.1) & (y_sub < y_cen + 0.1)
                bin_v1 = v1_sub[bin_mask]
                if len(bin_v1) > 0:
                    sums[sp_lbl][0][i] += np.sum(bin_v1)
                    sums[sp_lbl][1][i] += np.sum(bin_v1**2)
                    sums[sp_lbl][2][i] += len(bin_v1)
                    
    results = {}
    for sp_lbl in sums:
        s_v1, s_v1sq, count = sums[sp_lbl]
        y_centers = y_centers_dict[sp_lbl]
        mean_v1 = np.full_like(y_centers, np.nan)
        err_v1 = np.full_like(y_centers, np.nan)
        valid = count > 10
        if np.any(valid):
            mean_v1[valid] = s_v1[valid] / count[valid]
            var_v1 = np.maximum((s_v1sq[valid] / count[valid]) - (mean_v1[valid])**2, 0)
            err_v1[valid] = np.sqrt(var_v1) / np.sqrt(count[valid])
        results[sp_lbl] = (mean_v1, err_v1)
    return results

print("Processing ROOT files...")
res_lin = binned_mean_exact(ROOT_LIN, "Linear")
res_gau = binned_mean_exact(ROOT_GAU, "Gaussian")

# ── Plot ─────────────────────────────────────────────────────────────────────
print("Generating plot...")
fig, axs = plt.subplots(2, 3, figsize=(18, 12))

def get_slope(y_vals, v_vals):
    mask = (np.abs(y_vals) <= 0.8) & ~np.isnan(v_vals)
    if np.sum(mask) < 2: return np.nan
    sl, _, _, _, _ = linregress(y_vals[mask], v_vals[mask])
    return sl

def plot_panel(ax, sp_lbl, title, star_data, is_dv1=False, base_sp='p', ylim=None):
    sy, sv, serr = star_data
    if len(sy) > 0 and len(sv) > 0:
        ax.errorbar(sy, sv, yerr=serr, fmt='k*', mfc='black', mec='black', ecolor='black', label='STAR BES-I', markersize=12, zorder=10)
        s_slope = get_slope(sy, sv)
    else: s_slope = np.nan
    
    if is_dv1:
        v1p_l, e1p_l = res_lin[base_sp]
        v1m_l, e1m_l = res_lin[sp_lbl] # sp_lbl acts as anti-particle here
        v1_l = v1p_l - v1m_l
        err_l = np.sqrt(e1p_l**2 + e1m_l**2)
        
        v1p_g, e1p_g = res_gau[base_sp]
        v1m_g, e1m_g = res_gau[sp_lbl]
        v1_g = v1p_g - v1m_g
        err_g = np.sqrt(e1p_g**2 + e1m_g**2)
        
        y_cen = y_centers_dict[base_sp]
    else:
        v1_l, err_l = res_lin[sp_lbl]
        v1_g, err_g = res_gau[sp_lbl]
        y_cen = y_centers_dict[sp_lbl]
        
    valid_l = ~np.isnan(v1_l)
    ax.errorbar(y_cen[valid_l], v1_l[valid_l], yerr=err_l[valid_l], color='#2196F3', marker='s', ls='none', markersize=7, label='Linear', zorder=8)
    l_slope = get_slope(y_cen[valid_l], v1_l[valid_l])

    valid_g = ~np.isnan(v1_g)
    ax.errorbar(y_cen[valid_g], v1_g[valid_g], yerr=err_g[valid_g], color='#E91E63', marker='D', ls='none', markersize=7, label='Gaussian', zorder=9)
    g_slope = get_slope(y_cen[valid_g], v1_g[valid_g])
    
    text_str = ""
    if not np.isnan(s_slope): text_str += f"STAR slope: {s_slope:.3f}\n"
    if not np.isnan(g_slope): text_str += f"Gau. slope: {g_slope:.3f}\n"
    if not np.isnan(l_slope): text_str += f"Lin. slope: {l_slope:.3f}"
    
    if text_str:
        props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray')
        # Smart positioning based on expected v1 shape
        if 'Delta' in title or 'Pion' in title:
            loc = 'lower right'
        else:
            loc = 'upper left'
            
        if loc == 'upper left':
            ax.text(0.05, 0.95, text_str, transform=ax.transAxes, fontsize=11, va='top', bbox=props)
        else:
            ax.text(0.95, 0.05, text_str, transform=ax.transAxes, fontsize=11, ha='right', va='bottom', bbox=props)

    ax.axhline(0, color='gray', ls='--', lw=1.0)
    ax.axvline(0, color='gray', ls='--', lw=1.0)
    ax.set_xlabel(r'Rapidity $y$')
    ax.set_ylabel(r'Directed Flow $v_1$')
    ax.set_title(title)
    ax.set_xlim(-1.0, 1.0)
    if ylim:
        ax.set_ylim(ylim)
    ax.tick_params(direction='in', top=True, right=True)

# Row 1: Pions
plot_panel(axs[0, 0], 'pip', r'$\pi^+$ $v_1$', (s_y_pip, s_v1_pip, s_err_pip), ylim=[-0.04, 0.04])
plot_panel(axs[0, 1], 'pim', r'$\pi^-$ $v_1$', (s_y_pim, s_v1_pim, s_err_pim), ylim=[-0.04, 0.04])
plot_panel(axs[0, 2], 'pim', r'Pion $\Delta v_1$ ($\pi^+ - \pi^-$)', (y_centers_dict['pip'] if len(s_dv1_pi)>0 else [], s_dv1_pi, s_edv1_pi), is_dv1=True, base_sp='pip', ylim=[-0.04, 0.04])
axs[0, 0].legend(fontsize=10, loc='upper left')

# Row 2: Protons (Allow auto-scaling for pbar and delta v1 to show offset model points)
plot_panel(axs[1, 0], 'p', r'$p$ $v_1$', (s_y_p, s_v1_p, s_err_p), ylim=[-0.06, 0.06])
plot_panel(axs[1, 1], 'pbar', r'$\bar{p}$ $v_1$', (s_y_pbar, s_v1_pbar, s_err_pbar), ylim=None)
plot_panel(axs[1, 2], 'pbar', r'Proton $\Delta v_1$ ($p - \bar{p}$)', (y_centers_dict['p'] if len(s_dv1_p)>0 else [], s_dv1_p, s_edv1_p), is_dv1=True, base_sp='p', ylim=None)
axs[1, 0].legend(fontsize=10, loc='upper left')

fig.suptitle(r'Directed Flow $v_1(y)$ at $\sqrt{s_{NN}} = 7.7$ GeV (10-40%)', fontsize=18, fontweight='bold', y=1.02)
fig.tight_layout()
outpath = f'{OUTDIR}/exact_y_detailed_v1_all_models_v5.png'
fig.savefig(outpath, dpi=300, bbox_inches='tight')
print(f"\n✓ Saved: {outpath}")
