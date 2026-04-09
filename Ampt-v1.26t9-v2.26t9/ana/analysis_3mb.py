#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os

# -----------------------------
# User parameters
# -----------------------------
dat_file = "ampt.dat"  # path to AMPT output dat
output_dir = "plots"       # directory to save plots

# Particle PDG IDs
particle_dict = {
    'proton': 2212,
    'antiproton': -2212,
    'pion+': 211,
    'pion-': -211,
    'kaon+': 321,
    'kaon-': -321
}

# Particle-antiparticle pairs for combined plots
pairs = [
    ('proton', 'antiproton'),
    ('pion+', 'pion-'),
    ('kaon+', 'kaon-')
]

# -----------------------------
# Create output directory
# -----------------------------
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# -----------------------------
# Load AMPT data
# -----------------------------
data = []
with open(dat_file, 'r') as f:
    for line in f:
        if len(line.strip()) == 0:
            continue
        tokens = line.split()
        try:
            pdg = int(tokens[0])
            px, py, pz, E = map(float, tokens[1:5])
            data.append([pdg, px, py, pz, E])
        except ValueError:
            continue  # skip non-numeric lines

data = np.array(data)
print(f"Loaded {len(data)} particles from {dat_file}")

# -----------------------------
# Functions
# -----------------------------
def pT(px, py):
    return np.sqrt(px**2 + py**2)

def rapidity(E, pz):
    if E == pz:
        return 0.0
    return 0.5 * np.log((E + pz)/(E - pz))

def pseudorapidity(px, py, pz):
    p = np.sqrt(px**2 + py**2 + pz**2)
    if p - pz == 0:
        return 0.0
    return 0.5 * np.log((p + pz)/(p - pz))



# -----------------------------
# Individual particle analysis
# -----------------------------
particle_results = {}

for name, pdg_id in particle_dict.items():
    mask = data[:,0] == pdg_id
    particles = data[mask]
    
    if len(particles) == 0:
        print(f"No particles found for {name}")
        continue

    px = particles[:,1]
    py = particles[:,2]
    pz = particles[:,3]
    E  = particles[:,4]

    pt = pT(px, py)
    y  = np.array([rapidity(E[i], pz[i]) for i in range(len(E))])
    eta = np.array([pseudorapidity(px[i], py[i], pz[i]) for i in range(len(E))])

    particle_results[name] = {'pt': pt, 'y': y, 'eta': eta}

    # -------------------------
    # pT histogram
    plt.figure(figsize=(6,4))
    plt.hist(pt, bins=50, color='steelblue', alpha=0.8)
    plt.xlabel(r"$p_T$ [GeV/c]")
    plt.ylabel("Counts")
    plt.title(f"{name} pT distribution")
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{output_dir}/{name}_pT.png")
    plt.close()

    # -------------------------
    # Rapidity histogram
    plt.figure(figsize=(6,4))
    plt.hist(y, bins=50, color='seagreen', alpha=0.8)
    plt.xlabel("Rapidity y")
    plt.ylabel("Counts")
    plt.title(f"{name} Rapidity distribution")
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{output_dir}/{name}_rapidity.png")
    plt.close()

    # -------------------------
    # Pseudorapidity histogram
    plt.figure(figsize=(6,4))
    plt.hist(eta, bins=50, color='orange', alpha=0.8)
    plt.xlabel("Pseudorapidity η")
    plt.ylabel("Counts")
    plt.title(f"{name} Pseudorapidity distribution")
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{output_dir}/{name}_eta.png")
    plt.close()

    print(f"Individual plots saved for {name}: {len(particles)} particles")

# -----------------------------
# Combined particle-antiparticle plots
# -----------------------------
for p1, p2 in pairs:
    if p1 not in particle_results or p2 not in particle_results:
        print(f"Skipping combined plot for {p1} vs {p2} (data missing)")
        continue

    # pT comparison
    plt.figure(figsize=(6,4))
    plt.hist(particle_results[p1]['pt'], bins=50, alpha=0.6, label=p1)
    plt.hist(particle_results[p2]['pt'], bins=50, alpha=0.6, label=p2)
    plt.xlabel(r"$p_T$ [GeV/c]")
    plt.ylabel("Counts")
    plt.title(f"{p1} vs {p2} pT distribution")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{output_dir}/{p1}_vs_{p2}_pT.png")
    plt.close()

    # Rapidity comparison
    plt.figure(figsize=(6,4))
    plt.hist(particle_results[p1]['y'], bins=50, alpha=0.6, label=p1)
    plt.hist(particle_results[p2]['y'], bins=50, alpha=0.6, label=p2)
    plt.xlabel("Rapidity y")
    plt.ylabel("Counts")
    plt.title(f"{p1} vs {p2} Rapidity distribution")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{output_dir}/{p1}_vs_{p2}_rapidity.png")
    plt.close()

    # Pseudorapidity comparison
    plt.figure(figsize=(6,4))
    plt.hist(particle_results[p1]['eta'], bins=50, alpha=0.6, label=p1)
    plt.hist(particle_results[p2]['eta'], bins=50, alpha=0.6, label=p2)
    plt.xlabel("Pseudorapidity η")
    plt.ylabel("Counts")
    plt.title(f"{p1} vs {p2} Pseudorapidity distribution")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{output_dir}/{p1}_vs_{p2}_eta.png")
    plt.close()

    print(f"Combined plots saved for {p1} vs {p2}")

print("All analysis completed! Check the plots directory for results.")
