import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings("ignore")

# PDG Masses
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

def parse_ampt(filename):
    data = {pid: [] for pid in M0.keys()}
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
                            if abs(y) < 0.1:  # PRC exact cut
                                mt = np.sqrt(pt**2 + m**2)
                                data[pid].append(mt - m)
                except ValueError: pass
                finally: particles_left -= 1
                
    for k in data: data[k] = np.array(data[k])
    return data, max(num_events, 1)

files = {
    'Default': '../ana/ampt_default.dat',
    r'CQMF ($1\rho_0$)': '../ana/ampt_modified.dat',
    r'CQMF ($2\rho_0$)': '../ana/ampt_density2.dat',
    r'CQMF ($3\rho_0$)': '../ana/ampt_density3.dat'
}
colors = ['k', 'b', 'g', 'r']
markers = ['o', 's', '^', 'D']

datasets = {}
events = {}
for name, path in files.items():
    print(f"Loading {name}...")
    d, n = parse_ampt(path)
    datasets[name] = d
    events[name] = n

fig, axs = plt.subplots(1, 3, figsize=(15, 6), sharey=True)

species_groups = [
    ([211, -211], r'Pions ($\pi^\pm$)'),
    ([321, -321], r'Kaons ($K^\pm$)'),
    ([2212, -2212], r'Protons ($p, \bar{p}$)')
]

mt_bins = np.linspace(0.0, 2.5, 25)
mt_widths = mt_bins[1:] - mt_bins[:-1]
mt_centers = 0.5 * (mt_bins[1:] + mt_bins[:-1])
dy = 0.2  # y in [-0.1, 0.1]

for ax_idx, (pids, title) in enumerate(species_groups):
    ax = axs[ax_idx]
    
    for i, name in enumerate(files.keys()):
        d = datasets[name]
        n_ev = events[name]
        
        for pid in pids:
            h, _ = np.histogram(d[pid], bins=mt_bins)
            # Invariant yield: 1/(2*pi*mT) d2N/dydmT
            # mt_centers here is actually (mT-m0). 
            # So true mT = mt_centers + M0[pid]
            true_mt = mt_centers + M0[pid]
            
            with np.errstate(divide='ignore', invalid='ignore'):
                inv_yield = h / (n_ev * dy * mt_widths * 2 * np.pi * true_mt)
            
            # Offset negative particles for visibility (standard PRC technique)
            scale = 1.0
            if pid < 0: scale = 0.1
            
            label_name = f"{name} {LABELS[pid]}"
            if scale != 1.0: label_name += r" ($\times 10^{-1}$)"
            
            # Set open/closed markers for particles/antiparticles
            mfc = colors[i] if pid > 0 else 'none'
            ax.plot(mt_centers, inv_yield * scale, marker=markers[i], color=colors[i], 
                    markerfacecolor=mfc, linestyle='-', linewidth=1, markersize=6, label=label_name)

    ax.set_yscale('log')
    ax.set_xlabel(r'$m_T - m_0$ (GeV/$c^2$)')
    ax.set_title(title)
    if ax_idx == 0:
        ax.set_ylabel(r'$\frac{1}{2\pi m_T} \frac{d^2N}{dy dm_T}$ [$(GeV/c)^{-2}$]')
    ax.legend(fontsize=8, loc='upper right')
    # STAR BES limits
    ax.set_ylim(1e-5, 1e3)

fig.suptitle(r'STAR BES-I PRC 96 (2017) Figure 22 Analog: Invariant Yields at $\sqrt{s_{NN}} = 7.7$ GeV ($|y|<0.1$)', fontsize=14)
fig.tight_layout()
plt.savefig("fig22_mt_spectra.png", dpi=300)
print("Saved fig22_mt_spectra.png")
