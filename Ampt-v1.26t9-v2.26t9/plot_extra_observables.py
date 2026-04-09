import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
from scipy.stats import linregress
warnings.filterwarnings("ignore")

def extract_obs_data(filename):
    data = {}
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
                        mt = np.sqrt(pt**2 + m**2)
                        
                        if pid not in data:
                            data[pid] = {'pt':[], 'y':[], 'v1':[], 'v2':[], 'mt':[], 'm':m}
                        
                        data[pid]['pt'].append(pt)
                        data[pid]['y'].append(y)
                        data[pid]['v1'].append(v1)
                        data[pid]['v2'].append(v2)
                        data[pid]['mt'].append(mt)
                        
                except ValueError: pass
                finally: particles_left -= 1
    
    for pid in data:
        for k in ['pt', 'y', 'v1', 'v2', 'mt']:
            data[pid][k] = np.array(data[pid][k])
    return data

def calc_binned_mean(x_data, v_data, bins):
    v_mean = []
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    for i in range(len(bins)-1):
        mask = (x_data >= bins[i]) & (x_data < bins[i+1])
        if np.sum(mask) > 5: v_mean.append(np.mean(v_data[mask]))
        else: v_mean.append(np.nan)
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
    all_data[name] = extract_obs_data(fpath)

fig, axs = plt.subplots(2, 2, figsize=(15, 12))

# 1. Delta v2 (K+ - K-) Splitting
ax = axs[0, 0]
pt_bins = np.linspace(0.2, 2.0, 8)
for i, name in enumerate(files.keys()):
    if 321 in all_data[name] and -321 in all_data[name]:
        y_kp = all_data[name][321]['y']
        pt_kp = all_data[name][321]['pt']
        v2_kp_raw = all_data[name][321]['v2']
        mask_kp = (y_kp > -0.5) & (y_kp < 0.5)
        
        y_km = all_data[name][-321]['y']
        pt_km = all_data[name][-321]['pt']
        v2_km_raw = all_data[name][-321]['v2']
        mask_km = (y_km > -0.5) & (y_km < 0.5)
        
        _, v2_kp = calc_binned_mean(pt_kp[mask_kp], v2_kp_raw[mask_kp], pt_bins)
        _, v2_km = calc_binned_mean(pt_km[mask_km], v2_km_raw[mask_km], pt_bins)
        delta_v2 = v2_kp - v2_km
        
        ax.plot(pt_bins[:-1]+0.1, delta_v2, marker='o', color=colors[i], label=name)

ax.axhline(0, color='gray', linestyle='--')
ax.set_xlabel('$p_T$ (GeV/c)')
ax.set_ylabel('$\Delta v_2 (K^+ - K^-)$')
ax.set_title('Elliptic Flow Splitting vs $p_T$')
ax.legend()

# 2. dN/dy for Protons
ax = axs[0, 1]
y_bins = np.linspace(-3, 3, 20)
for i, name in enumerate(files.keys()):
    if 2212 in all_data[name]:
        y_data = all_data[name][2212]['y']
        hist, edges = np.histogram(y_data, bins=y_bins)
        bin_widths = edges[1:] - edges[:-1]
        # normalize by number of events (200) and bin width
        dndy = hist / (200.0 * bin_widths)
        ax.plot(edges[:-1]+(bin_widths/2), dndy, marker='s', color=colors[i], label=name)

ax.set_xlabel('Rapidity ($y$)')
ax.set_ylabel('$dN/dy$')
ax.set_title('Proton Rapidity Distribution')
ax.legend()

# 3. v1 slope dv1/dy at mid-rapidity for Protons
ax = axs[1, 0]
densities = [0, 1, 2, 3] # Default proxy is 0 
slopes = []
for name in files.keys():
    if 2212 in all_data[name]:
        y_data = all_data[name][2212]['y']
        v1_data = all_data[name][2212]['v1']
        # mid-rapidity mask
        mask = (y_data > -0.8) & (y_data < 0.8)
        if np.sum(mask) > 10:
            res = linregress(y_data[mask], v1_data[mask])
            slopes.append(res.slope)
        else:
            slopes.append(np.nan)

ax.plot(densities, slopes, marker='D', markersize=10, color='purple', linestyle='--', linewidth=2)
ax.set_xticks(densities)
ax.set_xticklabels(files.keys())
ax.set_ylabel('$dv_1/dy|_{y=0}$')
ax.set_title('Directed Flow Slope at Mid-Rapidity (Protons)')
ax.grid(True, alpha=0.3)

# 4. Transverse Mass spectra (mT - m0) for Kaons
ax = axs[1, 1]
mt_bins = np.linspace(0.0, 2.0, 20)
for i, name in enumerate(files.keys()):
    if 321 in all_data[name]:
        m0 = all_data[name][321]['m']
        mt_m0 = all_data[name][321]['mt'] - m0
        y_data = all_data[name][321]['y']
        mask = (y_data > -0.5) & (y_data < 0.5)
        
        hist, edges = np.histogram(mt_m0[mask], bins=mt_bins)
        bin_cen = 0.5 * (edges[:-1] + edges[1:])
        with np.errstate(divide='ignore', invalid='ignore'):
            inv_yield = hist / ((bin_cen + m0) * (edges[1]-edges[0]) * 200.0 * 1.0 * 2.0 * np.pi)
        ax.plot(bin_cen, inv_yield, marker='^', color=colors[i], label=name)

ax.set_yscale('log')
ax.set_xlabel('$m_T - m_0$ (GeV/c$^2$)')
ax.set_ylabel('$(1/2\pi m_T) d^2N/dydm_T$')
ax.set_title('$m_T$ Spectra offset for $K^+$')
ax.legend()

fig.suptitle('Advanced CQMF Observables in Au+Au @ 7.7 GeV', fontsize=18, y=0.98)
fig.tight_layout()
plt.savefig("advanced_observables.png", dpi=300)
print("Plot saved to advanced_observables.png")
