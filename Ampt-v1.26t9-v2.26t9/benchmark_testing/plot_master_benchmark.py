import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import json
import os
import csv

# ── Style ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size':   12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'legend.fontsize': 9,
    'lines.linewidth': 1.8,
    'lines.markersize': 5,
    'figure.facecolor': 'white',
})

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9'
DATA_FILE = os.path.join(BASE_DIR, 'benchmark_testing', 'model_v2_data.json')
CSV_FILE = os.path.join(BASE_DIR, 'benchmark_testing', 'benchmark_results.csv')
STAR_DIR = os.path.join(BASE_DIR, 'HEPData-ins1395151-v2-csv')
STAR_SPLIT_P = os.path.join(BASE_DIR, 'star_data', 'v2_splitting_p_pbar_7.7_10_40.csv')
OUTDIR = os.path.join(BASE_DIR, 'benchmark_testing')

# ── Load STAR Data ───────────────────────────────────────────────────────────
def load_star_csv(filename):
    filepath = os.path.join(STAR_DIR, filename)
    pt, v2, err = [], [], []
    with open(filepath) as f:
        for row in csv.reader(f):
            if not row or row[0].startswith('#') or row[0].startswith('PT') or row[0].startswith('$'): continue
            try:
                pt_val, v2_val = float(row[0]), row[1]
                if v2_val == '-': continue
                pt.append(pt_val)
                v2.append(float(v2_val))
                err.append(float(row[2]))
            except: continue
    return np.array(pt), np.array(v2), np.array(err)

def load_star_splitting(filepath):
    data = np.genfromtxt(filepath, delimiter=',', skip_header=1)
    return data[:, 0], data[:, 1], data[:, 4]

star_data = {
    'Pions': load_star_csv('Table110.csv'),
    'Kaons': load_star_csv('Table111.csv'),
    'Protons': load_star_csv('Table107.csv'),
    'Delta_v2': load_star_splitting(STAR_SPLIT_P)
}

# ── Model Configuration ──────────────────────────────────────────────────────
MODELS = {
    'M1_Default':    {'label': 'M1: Default AMPT',      'color': '#757575', 'marker': 'o', 'ls': '--'},
    'M2_Fixed_rho1': {'label': r'M2: Fixed $\rho_0$',   'color': '#FF9800', 'marker': 'v', 'ls': ':'},
    'M3_Fixed_rho2': {'label': r'M3: Fixed $2\rho_0$',  'color': '#F44336', 'marker': '^', 'ls': ':'},
    'M4_Fixed_rho3': {'label': r'M4: Fixed $3\rho_0$',  'color': '#9C27B0', 'marker': '<', 'ls': ':'},
    'M5_Linear':     {'label': 'M5: Linear Extrap.',    'color': '#2196F3', 'marker': 's', 'ls': '-.'},
    'M6_Gaussian':   {'label': 'M6: Gaussian Kernel',   'color': '#E91E63', 'marker': 'D', 'ls': '-'},
}

with open(DATA_FILE, 'r') as f:
    model_data = json.load(f)

# ── 1. Master Comparison Grid ────────────────────────────────────────────────
fig, axs = plt.subplots(2, 2, figsize=(15, 12))
axs = axs.flatten()

species_config = [
    ('Pions', r'Pion $v_2(p_T)$', r'$v_2$', [-0.02, 0.12]),
    ('Kaons', r'Kaon $v_2(p_T)$', r'$v_2$', [-0.02, 0.12]),
    ('Protons', r'Proton $v_2(p_T)$', r'$v_2$', [-0.02, 0.12]),
    ('Delta_v2', r'Proton-Antiproton Splitting $\Delta v_2$', r'$\Delta v_2$', [-0.05, 0.15])
]

for i, (spec, title, ylabel, ylim) in enumerate(species_config):
    ax = axs[i]
    
    # Plot Models
    for m_key, m_style in MODELS.items():
        if m_key not in model_data: continue
        d = model_data[m_key]
        
        if spec == 'Delta_v2':
            pt = np.array(d['Protons']['pt'])
            v2 = np.array(d['Protons']['v2']) - np.array(d['Antiprotons']['v2'])
            # Don't plot massive error bars for models here to keep it clean, or just plot lines
        else:
            pt = np.array(d[spec]['pt'])
            v2 = np.array(d[spec]['v2'])
            
        valid = ~np.isnan(v2) & (pt <= 2.2)
        ax.plot(pt[valid], v2[valid], color=m_style['color'], marker=m_style['marker'], 
                ls=m_style['ls'], lw=1.5, label=m_style['label'], zorder=5, alpha=0.8)

    # Plot STAR Data
    s_pt, s_v2, s_err = star_data[spec]
    ax.errorbar(s_pt, s_v2, yerr=s_err, fmt='k*', mfc='black', mec='black', ecolor='black',
                label='STAR BES-I (10-40%)', markersize=12, zorder=10)

    ax.axhline(0, color='black', ls='-', lw=0.5)
    ax.set_xlabel(r'$p_T$ (GeV/$c$)')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(0.0, 2.4)
    ax.set_ylim(ylim)
    if i == 0: ax.legend(fontsize=9, loc='upper left', ncol=2)
    ax.tick_params(direction='in', top=True, right=True)

fig.suptitle(r'Quantitative Benchmark: Model Comparisons vs STAR BES-I @ $\sqrt{s_{NN}} = 7.7$ GeV', fontsize=16, fontweight='bold', y=0.98)
fig.tight_layout()
out1 = os.path.join(OUTDIR, 'master_benchmark_grid.png')
fig.savefig(out1, dpi=300, bbox_inches='tight')

# ── 2. Global Ranking Bar Chart ──────────────────────────────────────────────
models, scores = [], []
with open(CSV_FILE, 'r') as f:
    reader = csv.reader(f)
    next(reader) # skip header
    for row in reader:
        models.append(row[0])
        scores.append(float(row[9]))

# Sort by score
sorted_indices = np.argsort(scores)
sorted_models = [models[i] for i in sorted_indices]
sorted_scores = [scores[i] for i in sorted_indices]

fig2, ax2 = plt.subplots(figsize=(10, 6))
# Define a color map for the models based on the MODELS dict
model_colors = [MODELS[m]['color'] for m in sorted_models]
bars = ax2.bar(sorted_models, sorted_scores, color=model_colors)

for bar in bars:
    yval = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2, yval + 2, round(yval, 1), ha='center', va='bottom', fontweight='bold')

ax2.set_ylabel(r'Global $\chi^2_{total}$ Score (Lower is Better)')
ax2.set_title(r'Simultaneous Benchmark Score across Pions, Kaons, Protons, and $\Delta v_2$')
ax2.grid(axis='y', linestyle='--', alpha=0.7)
ax2.set_xticklabels([MODELS[m]['label'].split(': ')[1] for m in sorted_models], rotation=25, ha='right')

fig2.tight_layout()
out2 = os.path.join(OUTDIR, 'global_ranking_score.png')
fig2.savefig(out2, dpi=300, bbox_inches='tight')

print(f"Generated {out1}")
print(f"Generated {out2}")
