import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings("ignore")

def extract_stopping_data(filename):
    data = {'p_y': [], 'pbar_y': [], 'kp_y': [], 'km_y': []}
    num_events = 0
    if not os.path.exists(filename): return data, num_events
    with open(filename, 'r') as f:
        particles_left = 0
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
                    px, py, pz, m = float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])
                    pt = np.sqrt(px**2 + py**2)
                    p_mag = np.sqrt(pt**2 + pz**2)
                    e = np.sqrt(p_mag**2 + m**2)
                    if pt > 0 and e > abs(pz):
                        y = 0.5 * np.log((e + pz) / (e - pz))
                        if pid == 2212: data['p_y'].append(y)
                        elif pid == -2212: data['pbar_y'].append(y)
                        elif pid == 321: data['kp_y'].append(y)
                        elif pid == -321: data['km_y'].append(y)
                except ValueError: pass
                finally: particles_left -= 1
    
    for k in data:
        data[k] = np.array(data[k])
    return data, max(num_events, 1)

files = {
    'Default': "ana/ampt_default.dat",
    'Mod ($\u03c1_0$)': "ana/ampt_modified.dat",
    'Mod (2$\u03c1_0$)': "ana/ampt_density2.dat",
    'Mod (3$\u03c1_0$)': "ana/ampt_density3.dat"
}
colors = ['royalblue', 'darkorange', 'forestgreen', 'firebrick']

all_data = {}
all_events = {}
for name, fpath in files.items():
    print(f"Loading {name}...")
    d, n_ev = extract_stopping_data(fpath)
    all_data[name] = d
    all_events[name] = n_ev

fig, axs = plt.subplots(1, 3, figsize=(18, 5))

y_bins = np.linspace(-3, 3, 20)
bin_cen = 0.5 * (y_bins[:-1] + y_bins[1:])
bin_widths = y_bins[1:] - y_bins[:-1]

for i, name in enumerate(files.keys()):
    d = all_data[name]
    n_ev = all_events[name]
    
    hist_p, _ = np.histogram(d['p_y'], bins=y_bins)
    hist_pbar, _ = np.histogram(d['pbar_y'], bins=y_bins)
    hist_kp, _ = np.histogram(d['kp_y'], bins=y_bins)
    hist_km, _ = np.histogram(d['km_y'], bins=y_bins)
    
    # 1. Net Proton dN/dy
    net_p = (hist_p - hist_pbar) / (n_ev * bin_widths)
    axs[0].plot(bin_cen, net_p, marker='s', color=colors[i], label=name)
    
    # 2. pbar/p ratio
    with np.errstate(divide='ignore', invalid='ignore'):
        pbar_p_ratio = hist_pbar / hist_p
    axs[1].plot(bin_cen, pbar_p_ratio, marker='o', color=colors[i], label=name)
    
    # 3. K+/K- ratio
    with np.errstate(divide='ignore', invalid='ignore'):
        kp_km_ratio = hist_kp / hist_km
    axs[2].plot(bin_cen, kp_km_ratio, marker='^', color=colors[i], label=name)

axs[0].axhline(0, color='gray', linestyle='--', linewidth=1)
axs[0].set_xlabel('Rapidity ($y$)')
axs[0].set_ylabel('Net-Proton $dN/dy$ ($p - \overline{p}$)')
axs[0].set_title('Baryon Stopping Profile')
axs[0].legend()

axs[1].axhline(1, color='gray', linestyle='--', linewidth=1)
axs[1].set_xlabel('Rapidity ($y$)')
axs[1].set_ylabel('$\overline{p} / p$ Ratio')
axs[1].set_title('Antimatter/Matter Asymmetry vs Rapidity')
axs[1].set_ylim(0, 1.5)
axs[1].legend()

axs[2].axhline(1, color='gray', linestyle='--', linewidth=1)
axs[2].set_xlabel('Rapidity ($y$)')
axs[2].set_ylabel('$K^+ / K^-$ Ratio')
axs[2].set_title('Strangeness Source Asymmetry vs Rapidity')
axs[2].legend()

fig.suptitle('Vector Potential Impacts: Baryon Stopping and Asymmetries', fontsize=15)
fig.tight_layout()
plt.savefig("baryon_stopping_asymmetry.png", dpi=300)
print("Plot saved to baryon_stopping_asymmetry.png")
