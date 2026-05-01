import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings("ignore")

M0 = {
    211: 0.13957, -211: 0.13957,
    321: 0.49367, -321: 0.49367,
    2212: 0.93827, -2212: 0.93827
}

LABELS = {
    211: r'$\pi^+$', -211: r'$\pi^-$',
    321: r'$K^+$', -321: r'$K^-$',
    2212: r'$p$', -2212: r'$\bar{p}$'
}

def parse_ampt_pt(filename):
    data = {pid: [] for pid in M0.keys()}
    if not os.path.exists(filename): return data
    
    particles_left = 0
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts: continue
            if particles_left == 0:
                try: particles_left = int(parts[2])
                except ValueError: pass
            else:
                try:
                    pid = int(parts[0])
                    if pid in M0:
                        px, py, pz = float(parts[1]), float(parts[2]), float(parts[3])
                        m = M0[pid]
                        pt = np.sqrt(px**2 + py**2)
                        p_mag = np.sqrt(pt**2 + pz**2)
                        e = np.sqrt(p_mag**2 + m**2)
                        if e > abs(pz):
                            y = 0.5 * np.log((e + pz)/(e - pz))
                            if abs(y) < 0.1:  # PRC exact cut
                                data[pid].append(pt)
                except ValueError: pass
                finally: particles_left -= 1
                
    for k in data: data[k] = np.array(data[k])
    return data

files = {
    'Default': '../ana/ampt_default.dat',
    r'CQMF ($1\rho_0$)': '../ana/ampt_modified.dat',
    r'CQMF ($2\rho_0$)': '../ana/ampt_density2.dat',
    r'CQMF ($3\rho_0$)': '../ana/ampt_density3.dat'
}

datasets = {}
for name, path in files.items():
    print(f"Loading {name}...")
    datasets[name] = parse_ampt_pt(path)

# Calculate mean pT and statistical errors
x_labels = list(files.keys())
x_pos = np.arange(len(x_labels))

mean_pt = {pid: [] for pid in M0.keys()}
err_pt = {pid: [] for pid in M0.keys()}

for name in x_labels:
    for pid in M0.keys():
        arr = datasets[name][pid]
        if len(arr) > 10:
            mean_pt[pid].append(np.mean(arr))
            err_pt[pid].append(np.std(arr)/np.sqrt(len(arr))) # Standard error of mean
        else:
            mean_pt[pid].append(np.nan)
            err_pt[pid].append(np.nan)

fig, axs = plt.subplots(1, 3, figsize=(15, 6), sharey=False)

groups = [
    ([211, -211], r'Pions', ['k', 'gray'], ['o', 's']),
    ([321, -321], r'Kaons', ['b', 'cyan'], ['o', 's']),
    ([2212, -2212], r'Protons', ['r', 'salmon'], ['o', 's'])
]

# Experimental 7.7 GeV STAR values (approx, from BES-I 0-5% centrality)
STAR_DATA = {
    211: 0.39, -211: 0.39,
    321: 0.61, -321: 0.60,
    2212: 0.79, -2212: 0.78
}

for ax_idx, (pids, title, colors, markers) in enumerate(groups):
    ax = axs[ax_idx]
    
    for i, pid in enumerate(pids):
        # Simulation line
        ax.errorbar(x_pos, mean_pt[pid], yerr=err_pt[pid], fmt=f'{markers[i]}-', 
                    color=colors[i], label=f'AMPT-CQMF {LABELS[pid]}', capsize=5, markersize=8)
        
        # STAR experiment band (horizontal line)
        if pid in STAR_DATA:
            ax.axhline(STAR_DATA[pid], color=colors[i], linestyle='--', alpha=0.6, 
                       label=f'STAR Exp {LABELS[pid]} (0-5%)')
    
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels)
    ax.set_ylabel(r'$\langle p_T \rangle$ (GeV/c)')
    ax.set_title(title)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

fig.suptitle(r'STAR BES-I PRC 96 (2017) Figure 18 Analog: Mean Transverse Momentum vs Vector Potential', fontsize=14)
fig.tight_layout()
plt.savefig("fig18_mean_pt.png", dpi=300)
print("Saved fig18_mean_pt.png")
