import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings("ignore")

def extract_flow_data(filename, target_pids):
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
    'Default': "ana_highstats/ampt_default.dat",
    'Mod ($\u03c1_0$)': "ana_highstats/ampt_modified.dat",
    'Mod (2$\u03c1_0$)': "ana_highstats/ampt_density2.dat",
    'Mod (3$\u03c1_0$)': "ana_highstats/ampt_density3.dat"
}
colors = ['royalblue', 'darkorange', 'forestgreen', 'firebrick']

all_p = {}
all_pi = {}
for name, fpath in files.items():
    print(f"Parsing {name} data...")
    all_p[name] = extract_flow_data(fpath, [2212])
    all_pi[name] = extract_flow_data(fpath, [211, -211])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

y_bins = np.linspace(-1.5, 1.5, 15)
pt_bins = np.linspace(0.2, 2.5, 12)

# Plot STAR data for comparison
def load_star_v2(filename, pt_col=0, v2_col=1, stat_err_col=4):
    if not os.path.exists(filename): return None, None, None
    df = pd.read_csv(filename)
    return df.iloc[:, pt_col].values, df.iloc[:, v2_col].values, df.iloc[:, stat_err_col].values

# Since we don't have pandas installed in the previous step, use numpy
def load_star_v2_np(filename, pt_col=0, v2_col=1, stat_err_col=4):
    if not os.path.exists(filename): return None, None, None
    try:
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        return data[:, pt_col], data[:, v2_col], data[:, stat_err_col]
    except: return None, None, None

# Plot v1 vs y (Protons)
for i, name in enumerate(files.keys()):
    if len(all_p[name]['y']) > 0:
        y_cen, v1, v1_err = calc_binned_flow(all_p[name]['y'], all_p[name]['v1'], y_bins)
        ax1.errorbar(y_cen, v1, yerr=v1_err, fmt='o-', color=colors[i], label=name, capsize=3, markersize=4)

# STAR v1 slope for protons at 7.7 GeV is around 0.015-0.02 (positive)
# We can add a representative slope line or just a note
ax1.set_ylim(-0.15, 0.15)
ax1.axhline(0, color='black', linewidth=0.8)
ax1.axvline(0, color='black', linewidth=0.8)
ax1.set_xlabel('Rapidity $y$', fontsize=12)
ax1.set_ylabel('Directed Flow $v_1$', fontsize=12)
ax1.set_title('Proton Directed Flow $v_1(y)$', fontsize=14)
ax1.legend(fontsize=10)

# Plot v2 vs pt (Pions)
for i, name in enumerate(files.keys()):
    if len(all_pi[name]['pt']) > 0:
        mask = (all_pi[name]['y'] > -0.5) & (all_pi[name]['y'] < 0.5)
        pt_cen, v2, v2_err = calc_binned_flow(all_pi[name]['pt'][mask], all_pi[name]['v2'][mask], pt_bins)
        ax2.errorbar(pt_cen, v2, yerr=v2_err, fmt='o-', color=colors[i], label=name, capsize=3, markersize=4)

# STAR data overlay for Pions
s_pt, s_v2, s_err = load_star_v2_np('star_data/v2_pip_7.7_0_80.csv')
if s_pt is not None:
    ax2.errorbar(s_pt, s_v2, yerr=s_err, fmt='ks', label='STAR (0-80%)', capsize=3, markersize=6, alpha=0.7)

ax2.set_xlabel('$p_T$ (GeV/c)', fontsize=12)
ax2.set_ylabel('Elliptic Flow $v_2$', fontsize=12)
ax2.set_title('Pion Elliptic Flow $v_2(p_T)$', fontsize=14)
ax2.legend(fontsize=10)

fig.suptitle('AMPT-CQMF vs STAR BES: Flow Observables (7.7 GeV Au+Au)', fontsize=16)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("publication_plots/flow_v1_v2_comparison.png", dpi=300)
print("Plot saved to publication_plots/flow_v1_v2_comparison.png")
