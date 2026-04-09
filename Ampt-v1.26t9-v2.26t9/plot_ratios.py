import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def parse_ampt_dat(filename):
    counts = {
        'pi+': 0, 'pi-': 0,
        'K+': 0, 'K-': 0,
        'p': 0, 'pbar': 0
    }
    pids = {
        211: 'pi+', -211: 'pi-',
        321: 'K+', -321: 'K-',
        2212: 'p', -2212: 'pbar'
    }
    
    if not os.path.exists(filename):
        return counts
        
    with open(filename, 'r') as f:
        particles_left = 0
        for line in f:
            parts = line.strip().split()
            if not parts: continue
            
            if particles_left == 0:
                try: particles_left = int(parts[2])
                except (ValueError, IndexError): pass
            else:
                try:
                    pid = int(parts[0])
                    if pid in pids: counts[pids[pid]] += 1
                except ValueError: pass
                finally: particles_left -= 1
    return counts

def calc_ratio(c1, c2):
    if c2 == 0: return 0, 0
    r = c1 / c2
    err = r * np.sqrt(1.0/c1 + 1.0/c2) if c1 > 0 else 0
    return r, err

files = {
    'Default': "ana/ampt_default.dat",
    'Mod ($\u03c1_0$)': "ana/ampt_modified.dat",
    'Mod (2$\u03c1_0$)': "ana/ampt_density2.dat",
    'Mod (3$\u03c1_0$)': "ana/ampt_density3.dat"
}

all_counts = {}
for name, fpath in files.items():
    print(f"Parsing {name} data...")
    all_counts[name] = parse_ampt_dat(fpath)

ratios = ['K+/pi+', 'K-/pi-', 'p/pi+', 'pbar/pi-']

x = np.arange(len(ratios))
width = 0.2

fig, ax = plt.subplots(figsize=(10, 6))

colors = ['royalblue', 'darkorange', 'forestgreen', 'firebrick']
offsets = [-1.5*width, -0.5*width, 0.5*width, 1.5*width]

for i, (name, counts) in enumerate(all_counts.items()):
    rs = [
        calc_ratio(counts['K+'], counts['pi+']),
        calc_ratio(counts['K-'], counts['pi-']),
        calc_ratio(counts['p'], counts['pi+']),
        calc_ratio(counts['pbar'], counts['pi-'])
    ]
    vals = [r[0] for r in rs]
    errs = [r[1] for r in rs]
    
    rects = ax.bar(x + offsets[i], vals, width, label=name, yerr=errs, capsize=3, color=colors[i])
    
    for rect, val in zip(rects, vals):
        height = rect.get_height()
        ax.annotate(f'{val:.3f}',
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 2),  
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=8, rotation=90)

ax.set_ylabel('Particle Ratio')
ax.set_title('Particle Ratios in Au+Au @ $\sqrt{s_{NN}}$=7.7 GeV vs Density')
ax.set_xticks(x)
ax.set_xticklabels(ratios)
ax.legend()
fig.tight_layout()
plt.savefig("particle_ratios_density_scan.png", dpi=300)
print("Plot saved to particle_ratios_density_scan.png")
