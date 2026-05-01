import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings("ignore")

# Figure 12 matches (pi+, pi-, K+, K-, p, pbar)
# Top row: pi+, K+, p
# Bottom row: pi-, K-, pbar

M0 = {
    211: 0.13957, -211: 0.13957,
    321: 0.49367, -321: 0.49367,
    2212: 0.93827, -2212: 0.93827
}

LAYOUT = {
    (0,0): (211, r'$\pi^+$', '(a)'),
    (1,0): (-211, r'$\pi^-$', '(b)'),
    (0,1): (321, r'$K^+$', '(c)'),
    (1,1): (-321, r'$K^-$', '(d)'),
    (0,2): (2212, r'$p$', '(e)'),
    (1,2): (-2212, r'$\bar{p}$', '(f)')
}

def parse_ampt_pt_spectra(filename):
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
                            if abs(y) < 0.1:  # Midrapidity cut to tightly match paper
                                data[pid].append(pt)
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

colors = ['r', 'b', 'g', 'm']
markers = ['o', 's', '^', 'D']

datasets = {}
events = {}
for name, path in files.items():
    print(f"Loading {name}...")
    d, n = parse_ampt_pt_spectra(path)
    datasets[name] = d
    events[name] = n

# 2x3 Grid exactly matching Fig 12
fig, axs = plt.subplots(2, 3, figsize=(14, 9), sharex=True, sharey=False)
plt.subplots_adjust(wspace=0.15, hspace=0.0)

# Paper uses pT range 0 to ~2 GeV/c
pt_bins = np.linspace(0.0, 2.2, 25)
pt_widths = pt_bins[1:] - pt_bins[:-1]
pt_centers = 0.5 * (pt_bins[1:] + pt_bins[:-1])
dy = 0.2  # |y| < 0.1 -> width 0.2

for row in range(2):
    for col in range(3):
        ax = axs[row, col]
        pid, title, letter = LAYOUT[(row, col)]
        
        # Determine baseline normalization scaling factor explicitly used in paper
        # Because we only simulate 0-5%, we DO NOT apply the arbitrary 1/2, 1/4 layout scalings
        # to separate curves because our curves are modifications, not separate centralities!
        # But we do keep the axis logs and markers exactly matching the paper aesthetic.
        
        for i, name in enumerate(files.keys()):
            d = datasets[name]
            n_ev = events[name]
            
            arr = d[pid]
            h, _ = np.histogram(arr, bins=pt_bins)
            
            # Invariant yield formula matching y-axis exactly
            # 1/(2*pi*pT) d2N/dydpT
            with np.errstate(divide='ignore', invalid='ignore'):
                inv_yield = h / (n_ev * dy * pt_widths * 2 * np.pi * pt_centers)
                inv_yield[inv_yield == 0] = np.nan  # Prevent vertical log-scale plunging for empty bins
                
            mfc = colors[i] if i == 0 else 'none'  # Matches paper style (first set filled, others open)
            ax.plot(pt_centers, inv_yield, marker=markers[i], color=colors[i], 
                    markerfacecolor=mfc, linestyle='-', linewidth=1, markersize=4, label=name)

        ax.set_yscale('log')
        
        # Limits matching paper (approx depending on particle)
        if abs(pid) == 211:
            ax.set_ylim(1e-4, 4e2)
        elif abs(pid) == 321:
            ax.set_ylim(1e-4, 5e1)
        elif abs(pid) == 2212:
            ax.set_ylim(1e-4, 5e1)
            
        # Specific structural text
        ax.text(0.1, 0.1, letter, transform=ax.transAxes, fontsize=12)
        ax.text(0.8, 0.9, title, transform=ax.transAxes, fontsize=18)
        
        if row == 0 and col == 0:
            ax.text(0.3, 0.9, r"Au+Au 7.7 GeV", transform=ax.transAxes, fontsize=12)
        
        # Ticks structure
        ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
        if row == 1:
            ax.set_xlabel(r'$p_T$ (GeV/c)', fontsize=14)
        if col == 0:
            ax.set_ylabel(r'$\frac{1}{2\pi p_T} \frac{d^2N}{dy dp_T} [(\mathrm{GeV}/c)^{-2}]$', fontsize=14)
            
        if row == 0 and col == 1:
            ax.legend(loc='lower left', frameon=False, fontsize=10)

fig.suptitle('Figure 12 Reproduction Layout: Invariant $p_T$ Spectra ($|y|<0.1$) vs Density Modification', fontsize=16)
plt.savefig("fig12_style_spectra.png", dpi=300, bbox_inches='tight')
print("Saved fig12_style_spectra.png")
