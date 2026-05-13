#!/usr/bin/env python3
"""Quick diagnostic: compute mean epsilon_2 and mean v2 for each h value."""
import os, numpy as np, uproot

ANA = "/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/gaussian_correction/ana"
H_VALS = [0.5, 0.8, 1.0, 1.2, 1.5, 2.0]

# STAR BES-I reference (pion v2, 10-40%, 7.7 GeV)
star_pt  = np.array([0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7])
star_v2  = np.array([0.010, 0.021, 0.038, 0.051, 0.062, 0.068, 0.071, 0.073])

print("="*70)
print(f"{'h (fm)':>8}  {'<eps2>':>10}  {'err_eps2':>10}  {'<v2> pions':>12}  {'chi2/ndf':>10}")
print("="*70)

results = []

for h in H_VALS:
    conf = f"h{h}"
    
    # 1. Eccentricity
    ecc_file = os.path.join(ANA, f"eccentricity_{conf}.dat")
    if os.path.exists(ecc_file):
        data = np.loadtxt(ecc_file)
        if len(data.shape) == 1:
            ecc_vals = np.array([data[1]])
        else:
            ecc_vals = data[:, 1]
        mean_ecc = np.mean(ecc_vals)
        err_ecc = np.std(ecc_vals) / np.sqrt(len(ecc_vals))
    else:
        mean_ecc = -1
        err_ecc = 0
    
    # 2. v2 (pions, |y|<1)
    root_file = os.path.join(ANA, f"ampt_{conf}.root")
    pt_bins = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8])
    v2_sums = np.zeros(len(pt_bins)-1)
    v2_sq   = np.zeros(len(pt_bins)-1)
    counts  = np.zeros(len(pt_bins)-1)
    
    for batch in uproot.iterate(f"{root_file}:ampt",
                                ["pid","px","py","pz","mass"],
                                step_size="100 MB", library="np"):
        pid  = batch["pid"]
        px   = batch["px"]
        py   = batch["py"]
        pz   = batch["pz"]
        mass = batch["mass"]
        
        e  = np.sqrt(px**2 + py**2 + pz**2 + mass**2)
        denom = e - pz
        denom[denom == 0] = 1e-9
        y  = 0.5 * np.log((e + pz) / denom)
        pt = np.sqrt(px**2 + py**2)
        phi = np.arctan2(py, px)
        v2 = np.cos(2*phi)
        
        mask = (np.abs(pid) == 211) & (np.abs(y) < 1.0)
        pt_m = pt[mask]
        v2_m = v2[mask]
        
        for j in range(len(pt_bins)-1):
            bm = (pt_m >= pt_bins[j]) & (pt_m < pt_bins[j+1])
            if np.sum(bm) > 0:
                v2_sums[j] += np.sum(v2_m[bm])
                v2_sq[j]   += np.sum(v2_m[bm]**2)
                counts[j]  += np.sum(bm)
    
    v2_mean = np.where(counts > 0, v2_sums / counts, 0)
    v2_err  = np.where(counts > 0, np.sqrt((v2_sq/counts - v2_mean**2)/counts), 0)
    pt_centers = (pt_bins[:-1] + pt_bins[1:]) / 2
    
    # Mean v2 in 0.5-1.5 GeV range
    mid_mask = (pt_centers >= 0.4) & (pt_centers <= 1.6)
    mean_v2 = np.mean(v2_mean[mid_mask]) if np.sum(mid_mask) > 0 else 0
    
    # Chi2 against STAR (interpolate AMPT to STAR pT points)
    ampt_v2_at_star = np.interp(star_pt, pt_centers, v2_mean)
    chi2 = np.sum((ampt_v2_at_star - star_v2)**2 / (0.005**2))  # assume ~0.005 combined error
    ndf = len(star_pt)
    
    print(f"{h:>8.1f}  {mean_ecc:>10.4f}  {err_ecc:>10.4f}  {mean_v2:>12.5f}  {chi2/ndf:>10.2f}")
    results.append((h, mean_ecc, err_ecc, mean_v2, v2_mean, v2_err, pt_centers, chi2/ndf))

print("="*70)

# Find the best h
best = min(results, key=lambda x: x[-1])
print(f"\n>>> BEST FIT: h = {best[0]} fm  (chi2/ndf = {best[-1]:.2f})")
print(f"    <epsilon_2> = {best[1]:.4f}")
print(f"    <v2> (0.4-1.6 GeV) = {best[3]:.5f}")
