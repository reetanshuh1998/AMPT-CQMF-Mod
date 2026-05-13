import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json
import os
import csv
from scipy.stats import linregress

# ── Style ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size':   12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'legend.fontsize': 9,
    'lines.linewidth': 1.8,
    'lines.markersize': 6,
    'figure.facecolor': 'white',
})

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9'
DATA_FILE = os.path.join(BASE_DIR, 'benchmark_testing', 'model_v1_data.json')
STAR_DIR = os.path.join(BASE_DIR, 'HEPData-ins1277069-v1-csv')
OUTDIR = os.path.join(BASE_DIR, 'benchmark_testing')
SLOPE_FILE = os.path.join(OUTDIR, 'v1_slopes.txt')

MODELS = {
    'M1_Default':    {'label': 'M1: Default AMPT',      'color': '#757575', 'marker': 'o', 'ls': '--'},
    'M2_Fixed_rho1': {'label': r'M2: Fixed $\rho_0$',   'color': '#FF9800', 'marker': 'v', 'ls': ':'},
    'M3_Fixed_rho2': {'label': r'M3: Fixed $2\rho_0$',  'color': '#F44336', 'marker': '^', 'ls': ':'},
    'M4_Fixed_rho3': {'label': r'M4: Fixed $3\rho_0$',  'color': '#9C27B0', 'marker': '<', 'ls': ':'},
    'M5_Linear':     {'label': 'M5: Linear Extrap.',    'color': '#2196F3', 'marker': 's', 'ls': '-.'},
    'M6_Gaussian':   {'label': 'M6: Gaussian Kernel',   'color': '#E91E63', 'marker': 'D', 'ls': '-'},
}

# ── Parse STAR Data ──────────────────────────────────────────────────────────
def load_star_v1(table_name, target_energy=7.7, particle_str='proton'):
    filepath = os.path.join(STAR_DIR, table_name)
    data = []
    in_target_block = False
    correct_particle = False
    
    if not os.path.exists(filepath):
        return np.array([]), np.array([]), np.array([])
        
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
                
            if row[0].startswith('#') or row[0].startswith('$'):
                continue
                
            if in_target_block and correct_particle:
                try:
                    data.append([float(row[0]), float(row[1]), float(row[2])])
                except ValueError:
                    continue
                    
    data = np.array(data)
    if len(data) == 0:
        return np.array([]), np.array([]), np.array([])
    return data[:, 0], data[:, 1], data[:, 2]

star_p_y, star_p_v1, star_p_err = load_star_v1('Table3.csv', 7.7, 'proton')
star_pi_y, star_pi_v1, star_pi_err = load_star_v1('Table4.csv', 7.7, '\\pi^{+}')

with open(DATA_FILE, 'r') as f:
    model_data = json.load(f)

# ── Calculate Slopes ─────────────────────────────────────────────────────────
slopes = {}
def get_slope(y_vals, v1_vals):
    mask = (np.abs(y_vals) <= 0.8) & ~np.isnan(v1_vals)
    if np.sum(mask) < 2: return np.nan
    slope, _, _, _, _ = linregress(y_vals[mask], v1_vals[mask])
    return slope

# ── Plotting ─────────────────────────────────────────────────────────────────
fig, axs = plt.subplots(1, 3, figsize=(18, 5))

species_config = [
    (axs[0], 'pip', r'Pion $\pi^+$ $v_1(y)$', [-0.04, 0.04], star_pi_y, star_pi_v1, star_pi_err),
    (axs[1], 'kp',  r'Kaon $K^+$ $v_1(y)$', [-0.04, 0.04], [], [], []),
    (axs[2], 'p',   r'Proton $p$ $v_1(y)$', [-0.08, 0.08], star_p_y, star_p_v1, star_p_err),
]

for m_key in MODELS.keys(): slopes[m_key] = {}

for ax, sp_key, title, ylim, sy, sv1, serr in species_config:
    for m_key, m_style in MODELS.items():
        if m_key not in model_data: continue
        d = model_data[m_key][sp_key]
        
        y = np.array(d['y'])
        v1 = np.array(d['v1'])
        
        valid = ~np.isnan(v1)
        ax.plot(y[valid], v1[valid], color=m_style['color'], marker=m_style['marker'], 
                ls=m_style['ls'], lw=1.5, label=m_style['label'], alpha=0.8)
                
        # Calculate slope
        slopes[m_key][sp_key] = get_slope(y, v1)

    # Plot STAR data
    if len(sy) > 0:
        ax.errorbar(sy, sv1, yerr=serr, fmt='k*', mfc='black', mec='black', ecolor='black',
                    label='STAR BES-I (10-40%)', markersize=12, zorder=10)

    ax.axhline(0, color='black', ls='-', lw=0.8)
    ax.axvline(0, color='black', ls='--', lw=0.5)
    ax.set_xlabel(r'Rapidity $y$')
    ax.set_ylabel(r'Directed Flow $v_1$')
    ax.set_title(title)
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(ylim)
    if ax == axs[0]:
        ax.legend(fontsize=9, loc='upper left')
    ax.tick_params(direction='in', top=True, right=True)

fig.suptitle(r'Directed Flow $v_1(y)$ Benchmark across Models at $\sqrt{s_{NN}} = 7.7$ GeV', fontsize=16, fontweight='bold', y=1.02)
fig.tight_layout()

out1 = os.path.join(OUTDIR, 'v1_benchmark_all_models.png')
fig.savefig(out1, dpi=300, bbox_inches='tight')

# Save Slopes
with open(SLOPE_FILE, 'w') as f:
    f.write(f"{'Model':<15} | {'Pion dv1/dy':<15} | {'Kaon dv1/dy':<15} | {'Proton dv1/dy':<15}\n")
    f.write("-" * 65 + "\n")
    for m_key in sorted(MODELS.keys()):
        p_sl = slopes[m_key].get('pip', np.nan)
        k_sl = slopes[m_key].get('kp', np.nan)
        pr_sl = slopes[m_key].get('p', np.nan)
        f.write(f"{m_key:<15} | {p_sl:>10.5f}      | {k_sl:>10.5f}      | {pr_sl:>10.5f}\n")

print(f"Saved {out1}")
print(f"Saved slopes to {SLOPE_FILE}")
