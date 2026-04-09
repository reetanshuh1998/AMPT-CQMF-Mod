import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings("ignore")

def extract_flow_data(filename):
    data = {'p': {'pt':[], 'y':[], 'v1':[], 'v2':[]},
            'pbar': {'pt':[], 'y':[], 'v1':[], 'v2':[]}}
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
                    px, py, pz, m = float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])
                    pt = np.sqrt(px**2 + py**2)
                    p_mag = np.sqrt(pt**2 + pz**2)
                    e = np.sqrt(p_mag**2 + m**2)
                    if pt > 0 and e > abs(pz):
                        y = 0.5 * np.log((e + pz) / (e - pz))
                        v1 = px / pt
                        v2 = (px**2 - py**2) / (pt**2)
                        
                        target_key = None
                        if pid == 2212: target_key = 'p'
                        elif pid == -2212: target_key = 'pbar'
                        
                        if target_key:
                            data[target_key]['pt'].append(pt)
                            data[target_key]['y'].append(y)
                            data[target_key]['v1'].append(v1)
                            data[target_key]['v2'].append(v2)
                except ValueError: pass
                finally: particles_left -= 1
                
    for key in data:
        for obs in data[key]:
            data[key][obs] = np.array(data[key][obs])
    return data

def calc_binned_flow(x_data, v_data, bins):
    v_mean = []
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    for i in range(len(bins)-1):
        mask = (x_data >= bins[i]) & (x_data < bins[i+1])
        if np.sum(mask) > 5:
            v_mean.append(np.mean(v_data[mask]))
        else:
            v_mean.append(np.nan)
    return bin_centers, np.array(v_mean)

files = {
    'Default': "ana/ampt_default.dat",
    'Mod ($\u03c1_0$)': "ana/ampt_modified.dat",
    'Mod (2$\u03c1_0$)': "ana/ampt_density2.dat",
    'Mod (3$\u03c1_0$)': "ana/ampt_density3.dat"
}
colors = ['royalblue', 'darkorange', 'forestgreen', 'firebrick']

all_data = {}
for name, fpath in files.items():
    print(f"Loading {name}...")
    all_data[name] = extract_flow_data(fpath)

fig, axs = plt.subplots(1, 2, figsize=(14, 6))

pt_bins = np.linspace(0.2, 2.0, 10)

# Proton - Antiproton v2 splitting
ax1 = axs[0]
for i, name in enumerate(files.keys()):
    d_p = all_data[name]['p']
    d_pbar = all_data[name]['pbar']
    
    if len(d_p['pt']) > 0 and len(d_pbar['pt']) > 0:
        mask_p = np.abs(d_p['y']) < 0.5
        mask_pb = np.abs(d_pbar['y']) < 0.5
        
        pt_cen, v2_p = calc_binned_flow(d_p['pt'][mask_p], d_p['v2'][mask_p], pt_bins)
        _, v2_pbar = calc_binned_flow(d_pbar['pt'][mask_pb], d_pbar['v2'][mask_pb], pt_bins)
        
        delta_v2 = v2_p - v2_pbar
        ax1.plot(pt_cen, delta_v2, marker='D', color=colors[i], label=name)

ax1.axhline(0, color='gray', linestyle='--')
ax1.set_xlabel('$p_T$ (GeV/c)')
ax1.set_ylabel('$\Delta v_2 (p - \overline{p})$')
ax1.set_title('Proton - Antiproton Elliptic Flow Splitting vs $p_T$')
ax1.legend()

# For a second plot, we just plot v2(p) to ensure it's not simply zero
ax2 = axs[1]
for i, name in enumerate(files.keys()):
    d_p = all_data[name]['p']
    if len(d_p['pt']) > 0:
        mask_p = np.abs(d_p['y']) < 0.5
        pt_cen, v2_p = calc_binned_flow(d_p['pt'][mask_p], d_p['v2'][mask_p], pt_bins)
        ax2.plot(pt_cen, v2_p, marker='o', color=colors[i], label=name)

ax2.axhline(0, color='gray', linestyle='--')
ax2.set_xlabel('$p_T$ (GeV/c)')
ax2.set_ylabel('$v_2$ (Protons)')
ax2.set_title('Proton Elliptic Flow ($v_2$) vs $p_T$')
ax2.legend()

fig.suptitle('Proton/Antiproton Vector Potential Splitting Signatures', fontsize=15)
fig.tight_layout()
plt.savefig("proton_v2_splitting.png", dpi=300)
print("Plot saved to proton_v2_splitting.png")
