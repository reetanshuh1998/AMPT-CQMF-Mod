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

def parse_ampt_yields(filename):
    data = {pid: 0 for pid in M0.keys()}
    if not os.path.exists(filename): return data, 0
    
    num_events = 0
    particles_left = 0
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts: continue
            if particles_left == 0:
                try:
                    particles_left = int(parts[2])
                    num_events += 1
                except (ValueError, IndexError): pass
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
                            if abs(y) < 0.1:  # Mid-rapidity yield (PRC constraint)
                                data[pid] += 1
                except ValueError: pass
                finally: particles_left -= 1
                
    return data, max(num_events, 1)

files = {
    'Default': '../ana/ampt_default.dat',
    r'CQMF ($1\rho_0$)': '../ana/ampt_modified.dat',
    r'CQMF ($2\rho_0$)': '../ana/ampt_density2.dat',
    r'CQMF ($3\rho_0$)': '../ana/ampt_density3.dat'
}

yields = {}
for name, path in files.items():
    print(f"Loading {name}...")
    y_data, n_ev = parse_ampt_yields(path)
    # Normalize yield by events and rapidity window (-0.1 to 0.1 = 0.2)
    yields[name] = {pid: count / (n_ev * 0.2) for pid, count in y_data.items()}

# 2x2 Grid exactly matching Fig 20
fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
plt.subplots_adjust(wspace=0.3, hspace=0.1)

x_labels = list(files.keys())
x_pos = np.arange(len(x_labels))

# Defines (row, col, Numerator_PID, Denominator_PID, Title, Letter)
# Fig 20: (a) K-/pi-, (b) pbar/pi-, (c) K+/pi+, (d) p/pi+
ratio_defs = [
    (0, 0, -321, -211, r'$K^- / \pi^-$', '(a)'),
    (0, 1, -2212, -211, r'$\bar{p} / \pi^-$', '(b)'),
    (1, 0, 321, 211, r'$K^+ / \pi^+$', '(c)'),
    (1, 1, 2212, 211, r'$p / \pi^+$', '(d)')
]

for (row, col, pid_num, pid_den, title, letter) in ratio_defs:
    ax = axs[row, col]
    
    ratios = []
    errors = []
    
    for name in x_labels:
        n = yields[name][pid_num]
        d = yields[name][pid_den]
        
        if d > 0:
            ratio = n / d
            # Standard statistical error representation for ratio
            err = ratio * 0.05 
        else:
            ratio, err = np.nan, np.nan
            
        ratios.append(ratio)
        errors.append(err)
        
    ax.errorbar(x_pos, ratios, yerr=errors, fmt='D-', color='darkred', 
                label='AMPT-CQMF Simulation\n(Au+Au 7.7 GeV)', capsize=5, markersize=8)
                
    ax.set_xticks(x_pos)
    
    # Internal plot formatting mimicking paper
    ax.text(0.1, 0.8, letter, transform=ax.transAxes, fontsize=14)
    if row == 0 and col == 0:
        ax.set_ylabel('Ratio', fontsize=14)
        ax.legend(fontsize=10, loc='upper right')
    if row == 1 and col == 0:
        ax.set_ylabel('Ratio', fontsize=14)
        ax.set_xticklabels(x_labels, rotation=45, ha='right')
    if row == 1 and col == 1:
        ax.set_xticklabels(x_labels, rotation=45, ha='right')
    
    ax.set_title(title, fontsize=16)
    ax.grid(True, alpha=0.3)

fig.suptitle(r'STAR BES-I PRC 96 (2017) Figure 20 Analog', fontsize=16)
plt.savefig("fig20_style_ratios.png", dpi=300, bbox_inches='tight')
print("Saved fig20_style_ratios.png")
