import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings("ignore")

def parse_ampt(filename):
    """Extract proton, antiproton, K+, K- four-momenta from ampt.dat"""
    data = {
        2212:  {'pt': [], 'y': []},   # proton
        -2212: {'pt': [], 'y': []},   # antiproton
        321:   {'pt': [], 'y': []},   # K+
        -321:  {'pt': [], 'y': []},   # K-
        211:   {'pt': [], 'y': []},   # pi+
        -211:  {'pt': [], 'y': []},   # pi-
    }
    if not os.path.exists(filename):
        return data, 0
    
    num_events = 0
    particles_left = 0
    
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            if particles_left == 0:
                try:
                    particles_left = int(parts[2])
                    num_events += 1
                except (ValueError, IndexError):
                    pass
            else:
                try:
                    pid = int(parts[0])
                    if pid in data:
                        px = float(parts[1])
                        py = float(parts[2])
                        pz = float(parts[3])
                        m  = float(parts[4])
                        pt = np.sqrt(px**2 + py**2)
                        e  = np.sqrt(pt**2 + pz**2 + m**2)
                        if e > abs(pz) and pt > 0:
                            y = 0.5 * np.log((e + pz) / (e - pz))
                            data[pid]['pt'].append(pt)
                            data[pid]['y'].append(y)
                except (ValueError, IndexError):
                    pass
                finally:
                    particles_left = max(0, particles_left - 1)
    
    for pid in data:
        for k in data[pid]:
            data[pid][k] = np.array(data[pid][k])
    return data, max(num_events, 1)

# ---- Configuration ----
files = {
    'Default':               'ana/ampt_default.dat',
    r'Mod ($\rho_0$)':       'ana/ampt_modified.dat',
    r'Mod ($2\rho_0$)':      'ana/ampt_density2.dat',
    r'Mod ($3\rho_0$)':      'ana/ampt_density3.dat',
}
colors = ['royalblue', 'darkorange', 'forestgreen', 'firebrick']

all_data = {}
all_nevt = {}
for label, fp in files.items():
    print(f"Loading {label}...")
    d, n = parse_ampt(fp)
    all_data[label] = d
    all_nevt[label] = n

# ---- Plot ----
fig, axs = plt.subplots(2, 6, figsize=(30, 10))

pt_bins = np.linspace(0.0, 3.0, 25)
pt_cen  = 0.5 * (pt_bins[1:] + pt_bins[:-1])
pt_w    = pt_bins[1:] - pt_bins[:-1]

y_bins  = np.linspace(-3.0, 3.0, 25)
y_cen   = 0.5 * (y_bins[1:] + y_bins[:-1])
y_w     = y_bins[1:] - y_bins[:-1]

species = [
    (2212,  'Proton ($p$)'),
    (-2212, r'Antiproton ($\bar{p}$)'),
    (321,   'Kaon ($K^+$)'),
    (-321,  'Kaon ($K^-$)'),
    (211,   r'Pion ($\pi^+$)'),
    (-211,  r'Pion ($\pi^-$)'),
]

# Row 0: pT spectra   |   Row 1: rapidity distributions
for col, (pid, plabel) in enumerate(species):
    ax_pt = axs[0, col]
    ax_y  = axs[1, col]
    
    for i, label in enumerate(files.keys()):
        n_ev = all_nevt[label]
        pt_arr = all_data[label][pid]['pt']
        y_arr  = all_data[label][pid]['y']
        
        # pT spectrum (dN/dpT, mid-rapidity |y|<0.5)
        mid_mask = np.abs(y_arr) < 0.5
        h_pt, _ = np.histogram(pt_arr[mid_mask], bins=pt_bins)
        dN_dpt  = h_pt / (n_ev * pt_w)
        ax_pt.plot(pt_cen, dN_dpt, marker='o', markersize=4, color=colors[i], label=label)
        
        # Rapidity distribution (dN/dy, all pT)
        h_y, _ = np.histogram(y_arr, bins=y_bins)
        dN_dy  = h_y / (n_ev * y_w)
        ax_y.plot(y_cen, dN_dy, marker='s', markersize=4, color=colors[i], label=label)
    
    ax_pt.set_xlabel('$p_T$ (GeV/c)')
    ax_pt.set_ylabel('$dN/dp_T$ ($|y|<0.5$)')
    ax_pt.set_title(f'{plabel} — $p_T$ Spectrum')
    ax_pt.set_yscale('log')
    ax_pt.legend(fontsize=7)
    
    ax_y.set_xlabel('Rapidity ($y$)')
    ax_y.set_ylabel('$dN/dy$')
    ax_y.set_title(f'{plabel} — Rapidity Distribution')
    ax_y.legend(fontsize=7)

fig.suptitle(r'Particle Production: $p_T$ Spectra and Rapidity Distributions vs Density',
             fontsize=16, fontweight='bold')
fig.tight_layout()
plt.savefig('proton_kaon_production.png', dpi=300)
print("Saved: proton_kaon_production.png")
