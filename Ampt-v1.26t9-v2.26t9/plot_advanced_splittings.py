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
    'Default (No Medium)':    "local_density_approach/ana/ampt_default.dat",
    'Fixed ρ=1ρ₀':           "local_density_approach/ana/ampt_fixed_rho1.dat",
    'Fixed ρ=2ρ₀':           "local_density_approach/ana/ampt_fixed_rho2.dat",
    'Fixed ρ=3ρ₀':           "local_density_approach/ana/ampt_fixed_rho3.dat",
    'Local Density (iqmc=2)': "local_density_approach/ana/ampt_localdensity_fixed.dat",
}
colors = ['royalblue', 'darkorange', 'forestgreen', 'firebrick', 'purple']

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

# STAR data overlay for comparison
def load_star_data_np(filename, x_col=0, y_col=1, err_col=4):
    if not os.path.exists(filename): return None, None, None
    try:
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        return data[:, x_col], data[:, y_col], data[:, err_col]
    except: return None, None, None

# Plot STAR v2 splitting
s_pt, s_dv2, s_err = load_star_data_np('star_data/v2_splitting_p_pbar_7.7_10_40.csv')
if s_pt is not None:
    ax1.errorbar(s_pt, s_dv2, yerr=s_err, fmt='ks', label='STAR (10-40%)', capsize=3, markersize=6, alpha=0.8)

ax1.axhline(0, color='black', linewidth=0.8)
ax1.set_xlabel('$p_T$ (GeV/c)', fontsize=12)
ax1.set_ylabel(r'$\Delta v_2 (p - \overline{p})$', fontsize=12)
ax1.set_title('Proton - Antiproton Elliptic Flow Splitting', fontsize=14)
ax1.legend(fontsize=10)

# Plot simulated v2(p) for all models on ax2
ax2 = axs[1]
for i, name in enumerate(files.keys()):
    d_p = all_data[name]['p']
    if len(d_p['pt']) > 0:
        mask_p = np.abs(d_p['y']) < 0.5
        pt_cen, v2_p = calc_binned_flow(d_p['pt'][mask_p], d_p['v2'][mask_p], pt_bins)
        ax2.plot(pt_cen, v2_p, marker='D', color=colors[i], label=name)

# Plot STAR v2(p) for comparison on ax2
s_pt_p, s_v2_p, s_err_p = load_star_data_np('star_data/v2_proton_7.7_0_10.csv')
if s_pt_p is not None:
    ax2.errorbar(s_pt_p, s_v2_p, yerr=s_err_p, fmt='ks', label='STAR (0-10%)', capsize=3, markersize=6, alpha=0.8)
    ax2.errorbar(s_pt_p, s_v2_p, yerr=s_err_p, fmt='ks', label='STAR (0-10%)', capsize=3, markersize=6, alpha=0.8)

ax2.axhline(0, color='black', linewidth=0.8)
ax2.set_xlabel('$p_T$ (GeV/c)', fontsize=12)
ax2.set_ylabel('$v_2$ (Protons)', fontsize=12)
ax2.set_title('Proton Elliptic Flow $v_2(p_T)$', fontsize=14)
ax2.legend(fontsize=10)

fig.suptitle('Proton/Antiproton Vector Potential Signatures (7.7 GeV Au+Au)', fontsize=16)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("publication_plots/proton_v2_splitting_comparison.png", dpi=300)
print("Plot saved to publication_plots/proton_v2_splitting_comparison.png")
