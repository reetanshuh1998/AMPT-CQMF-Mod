import numpy as np
import matplotlib.pyplot as plt
import os

def extract_pt_data(filename, target_pids):
    pt_list = []
    if not os.path.exists(filename): return pt_list
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
                        px, py = float(parts[1]), float(parts[2])
                        pt = np.sqrt(px**2 + py**2)
                        pt_list.append(pt)
                except ValueError: pass
                finally: particles_left -= 1
    return np.array(pt_list)

files = {
    'Default': "ana/ampt_default.dat",
    'Mod ($\u03c1_0$)': "ana/ampt_modified.dat",
    'Mod (2$\u03c1_0$)': "ana/ampt_density2.dat",
    'Mod (3$\u03c1_0$)': "ana/ampt_density3.dat"
}
colors = ['royalblue', 'darkorange', 'forestgreen', 'firebrick']

all_pt_k = {}
all_pt_pi = {}
for name, fpath in files.items():
    print(f"Parsing {name} data...")
    all_pt_k[name] = extract_pt_data(fpath, [321, -321])
    all_pt_pi[name] = extract_pt_data(fpath, [211, -211])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
bins = np.linspace(0, 3.0, 30)

bin_centers = 0.5*(bins[1:] + bins[:-1])

for i, name in enumerate(files.keys()):
    h_k, _ = np.histogram(all_pt_k[name], bins=bins, density=True)
    ax1.plot(bin_centers, h_k, color=colors[i], label=name, marker='o', markersize=4)
    
    h_pi, _ = np.histogram(all_pt_pi[name], bins=bins, density=True)
    ax2.plot(bin_centers, h_pi, color=colors[i], label=name, marker='s', markersize=4)

ax1.set_yscale('log')
ax1.set_xlabel('$p_T$ (GeV/c)')
ax1.set_ylabel('$(1/N) dN/dp_T$')
ax1.set_title('Normalized $p_T$ Spectra: Kaons ($K^\pm$)')
ax1.legend()

ax2.set_yscale('log')
ax2.set_xlabel('$p_T$ (GeV/c)')
ax2.set_ylabel('$(1/N) dN/dp_T$')
ax2.set_title('Normalized $p_T$ Spectra: Pions ($\pi^\pm$)')
ax2.legend()

fig.suptitle('Transverse Momentum ($p_T$) Spectra vs Density', fontsize=14)
fig.tight_layout()
plt.savefig("pt_spectra_density_scan.png", dpi=300)
print("Plot saved to pt_spectra_density_scan.png")
