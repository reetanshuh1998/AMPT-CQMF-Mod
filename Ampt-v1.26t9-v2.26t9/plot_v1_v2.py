import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings("ignore")

def extract_flow_data(filename, target_pids, y_cut=None):
    """Extract flow variables.  Apply |y| < y_cut if y_cut is given."""
    data = {'pt': [], 'y': [], 'v1': [], 'v2': []}
    if not os.path.exists(filename): return data
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
                    if pid in target_pids:
                        px, py, pz, m = float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])
                        pt = np.sqrt(px**2 + py**2)
                        p_mag = np.sqrt(pt**2 + pz**2)
                        e = np.sqrt(p_mag**2 + m**2)
                        if pt > 0 and e > abs(pz):
                            y = 0.5 * np.log((e + pz) / (e - pz))
                            if y_cut is not None and abs(y) >= y_cut:
                                pass  # outside acceptance
                            else:
                                v1 = px / pt
                                v2 = (px**2 - py**2) / (pt**2)
                                data['pt'].append(pt)
                                data['y'].append(y)
                                data['v1'].append(v1)
                                data['v2'].append(v2)
                except ValueError: pass
                finally: particles_left -= 1
    for k in data: data[k] = np.array(data[k])
    return data

def calc_binned_flow(x_data, v_data, bins):
    v_mean, v_err = [], []
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    for i in range(len(bins)-1):
        mask = (x_data >= bins[i]) & (x_data < bins[i+1])
        if np.sum(mask) > 10:
            v_mean.append(np.mean(v_data[mask]))
            v_err.append(np.std(v_data[mask]) / np.sqrt(np.sum(mask)))
        else:
            v_mean.append(np.nan)
            v_err.append(np.nan)
    return bin_centers, np.array(v_mean), np.array(v_err)

files = {
    'Default': "ana/ampt_default.dat",
    'Mod ($\u03c1_0$)': "ana/ampt_modified.dat",
    'Mod (2$\u03c1_0$)': "ana/ampt_density2.dat",
    'Mod (3$\u03c1_0$)': "ana/ampt_density3.dat"
}
colors = ['royalblue', 'darkorange', 'forestgreen', 'firebrick']

all_p = {}
all_pi = {}
for name, fpath in files.items():
    print(f"Parsing {name} data...")
    # Protons: full rapidity range for v1(y) plot
    all_p[name] = extract_flow_data(fpath, [2212])
    # Pions: restrict to |y| < 1.0 for mid-rapidity v2(pT) plot
    all_pi[name] = extract_flow_data(fpath, [211, -211], y_cut=1.0)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

y_bins = np.linspace(-1.5, 1.5, 15)
pt_bins = np.linspace(0.2, 2.5, 12)

for i, name in enumerate(files.keys()):
    if len(all_p[name]['y']) > 0:
        y_cen, v1, v1_err = calc_binned_flow(all_p[name]['y'], all_p[name]['v1'], y_bins)
        ax1.errorbar(y_cen, v1, yerr=v1_err, fmt='o-', color=colors[i], label=name, capsize=3)
    
    if len(all_pi[name]['pt']) > 0:
        pt_cen, v2, v2_err = calc_binned_flow(all_pi[name]['pt'], all_pi[name]['v2'], pt_bins)
        ax2.errorbar(pt_cen, v2, yerr=v2_err, fmt='o-', color=colors[i], label=name, capsize=3)

ax1.axhline(0, color='gray', linestyle='--', linewidth=1)
ax1.set_xlabel('Rapidity ($y$)')
ax1.set_ylabel('Directed Flow $v_1$')
ax1.set_title('$v_1$ vs Rapidity (Protons)')
ax1.legend()

ax2.axhline(0, color='gray', linestyle='--', linewidth=1)
ax2.set_xlabel('$p_T$ (GeV/c)')
ax2.set_ylabel('Elliptic Flow $v_2$')
ax2.set_title(r'$v_2$ vs $p_T$ (Pions, $|y|<1.0$)')
ax2.legend()

fig.suptitle('Directed Flow ($v_1$) and Elliptic Flow ($v_2$) vs Density', fontsize=15)
fig.tight_layout()
plt.savefig("flow_v1_v2_density_scan.png", dpi=300)
print("Plot saved to flow_v1_v2_density_scan.png")
