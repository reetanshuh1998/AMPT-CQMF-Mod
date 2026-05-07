"""
Publication-quality dN/dy + Ratio plots — WITH PROPER KINEMATIC CUTS
Au+Au @ 7.7 GeV — 200k events

The sharp proton peaks near y_beam = ±2.1 are SPECTATOR NUCLEONS.
They carry no medium physics. Real experiments (STAR) cut them out using:
  - pT > 0.4 GeV/c  (removes cold spectators)
  - |y| < 1.0       (mid-rapidity acceptance)

This script applies the same cuts.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, time
import uproot

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 11,
    'axes.labelsize': 13,
    'axes.titlesize': 14,
    'legend.fontsize': 8.5,
    'figure.dpi': 200,
    'lines.linewidth': 1.8,
})

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'local_density_approach', 'ana')
OUTDIR = 'publication_plots'
os.makedirs(OUTDIR, exist_ok=True)

CONFIGS = {
    'Default (No Medium)': os.path.join(BASE, 'ampt_default.root'),
    r'Fixed $\rho=1\rho_0$': os.path.join(BASE, 'ampt_fixed_rho1.root'),
    r'Fixed $\rho=2\rho_0$': os.path.join(BASE, 'ampt_fixed_rho2.root'),
    r'Fixed $\rho=3\rho_0$': os.path.join(BASE, 'ampt_fixed_rho3.root'),
    'Local Density Model': os.path.join(BASE, 'ampt_localdensity.root'),
}
CONFIG_NAMES = list(CONFIGS.keys())

COLORS = {
    'Default (No Medium)': 'royalblue',
    r'Fixed $\rho=1\rho_0$': 'darkorange',
    r'Fixed $\rho=2\rho_0$': 'forestgreen',
    r'Fixed $\rho=3\rho_0$': 'firebrick',
    'Local Density Model': 'purple'
}
MARKERS = {
    'Default (No Medium)': 'o',
    r'Fixed $\rho=1\rho_0$': 's',
    r'Fixed $\rho=2\rho_0$': '^',
    r'Fixed $\rho=3\rho_0$': 'D',
    'Local Density Model': 'P'
}
LS = {
    'Default (No Medium)': '-',
    r'Fixed $\rho=1\rho_0$': '--',
    r'Fixed $\rho=2\rho_0$': '-.',
    r'Fixed $\rho=3\rho_0$': ':',
    'Local Density Model': (0, (3, 1, 1, 1))
}

PID = dict(pip=211, pim=-211, kp=321, km=-321, pr=2212, pbar=-2212)

# ── Parse ────────────────────────────────────────────────────────────────
def parse_root(filepath):
    if not os.path.exists(filepath):
        return 0, None, None, None
    t0 = time.time()
    with uproot.open(filepath) as f:
        tree = f["ampt"]
        data = tree.arrays(["b", "pid", "px", "py", "pz", "mass"], library="np")
    b = data["b"]
    nevents = len(np.unique(b))
    pid = data["pid"]
    px, py, pz, mass = data["px"], data["py"], data["pz"], data["mass"]

    pt = np.sqrt(px**2 + py**2)
    e = np.sqrt(px**2 + py**2 + pz**2 + mass**2)
    denom = e - pz
    ok = (denom > 1e-9) & (e > 1e-9)
    rap = np.full_like(e, np.nan)
    rap[ok] = 0.5 * np.log((e[ok] + pz[ok]) / denom[ok])

    print(f"  {os.path.basename(filepath)}: {nevents} ev, {len(pid)} particles ({time.time()-t0:.1f}s)")
    return nevents, pid, rap, pt

print("Parsing all ROOT files...")
DATA = {}
for name, fpath in CONFIGS.items():
    print(f"[{name}]")
    nev, pid, rap, pt = parse_root(fpath)
    if nev > 0:
        DATA[name] = {'nev': nev, 'pid': pid, 'y': rap, 'pt': pt}

y_beam = np.arccosh(7.7 / (2 * 0.938))

species_list = [
    (r'Pions ($\pi^{\pm}$)',      [PID['pip'], PID['pim']]),
    (r'Kaons ($K^{\pm}$)',        [PID['kp'],  PID['km']]),
    (r'Protons ($p + \bar{p}$)',  [PID['pr'],  PID['pbar']]),
]


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 1: dN/dy — WITH pT CUTS
# ═══════════════════════════════════════════════════════════════════════
print("\nFigure 1: dN/dy (with pT-cuts)...")

ybins = np.arange(-4.0, 4.01, 0.25)
bc = 0.5 * (ybins[:-1] + ybins[1:])
bw = ybins[1] - ybins[0]

fig1, axs1 = plt.subplots(1, 3, figsize=(18, 5.5))

for idx, (sp_label, sp_pids) in enumerate(species_list):
    ax = axs1[idx]
    pt_cut = 0.4 if 'Proton' in sp_label else 0.2
    
    for name in CONFIG_NAMES:
        if name not in DATA:
            continue
        D = DATA[name]
        nev = D['nev']
        mask = np.isin(D['pid'], sp_pids) & (D['pt'] > pt_cut)
        h, _ = np.histogram(D['y'][mask], bins=ybins)
        dNdy = h / (nev * bw)
        err = np.sqrt(h) / (nev * bw)
        ax.errorbar(bc, dNdy, yerr=err,
                    color=COLORS[name], marker=MARKERS[name], ls=LS[name],
                    markersize=4, label=name, capsize=2, elinewidth=0.7)

    ax.set_xlabel('Rapidity $y$')
    ax.set_ylabel(r'$(1/N_{\mathrm{ev}})\, dN/dy$')
    ax.set_title(sp_label + f' ($p_T > {pt_cut}$ GeV/$c$)')
    ax.set_xlim(-4, 4)
    ax.set_ylim(bottom=0)
    ax.axvline(0, color='gray', lw=0.5, ls='--', alpha=0.4)
    ax.legend(fontsize=7, loc='upper right')

    if 'Proton' in sp_label:
        ax.axvline(y_beam, color='red', lw=0.7, ls=':', alpha=0.4)
        ax.axvline(-y_beam, color='red', lw=0.7, ls=':', alpha=0.4)

fig1.suptitle(r'Au+Au @ 7.7 GeV — Rapidity Distributions (with $p_T$ cuts)',
              fontsize=15, fontweight='bold')
fig1.tight_layout()
fig1.savefig(os.path.join(OUTDIR, 'fig1_dndy.png'), dpi=200, bbox_inches='tight')
print("  Saved fig1_dndy.png")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 2: Baryon stopping (with pT > 0.4 cut)
# ═══════════════════════════════════════════════════════════════════════
print("Figure 2: Baryon stopping (pT > 0.4)...")
fig2, axs2 = plt.subplots(1, 2, figsize=(14, 5.5))

for name in CONFIG_NAMES:
    if name not in DATA:
        continue
    D = DATA[name]
    nev, pid, y, pt = D['nev'], D['pid'], D['y'], D['pt']
    kw = dict(color=COLORS[name], marker=MARKERS[name], ls=LS[name],
              markersize=3.5, capsize=1.5, elinewidth=0.6)

    # Apply pT > 0.4 cut
    mask_p   = (pid == PID['pr'])   & (pt > 0.4)
    mask_pb  = (pid == PID['pbar']) & (pt > 0.4)
    mask_kp  = (pid == PID['kp'])   & (pt > 0.2)
    mask_km  = (pid == PID['km'])   & (pt > 0.2)

    hp, _  = np.histogram(y[mask_p],  bins=ybins)
    hpb, _ = np.histogram(y[mask_pb], bins=ybins)
    hkp, _ = np.histogram(y[mask_kp], bins=ybins)
    hkm, _ = np.histogram(y[mask_km], bins=ybins)

    net_p = (hp - hpb) / (nev * bw)
    err_np = np.sqrt(hp + hpb) / (nev * bw)
    axs2[0].errorbar(bc, net_p, yerr=err_np, label=name, **kw)

    with np.errstate(divide='ignore', invalid='ignore'):
        krat = np.where(hkm > 0, hkp / hkm, np.nan)
        krat_err = np.where(hkm > 0,
                            krat * np.sqrt(1.0/np.maximum(hkp, 1) + 1.0/np.maximum(hkm, 1)),
                            np.nan)
    axs2[1].errorbar(bc, krat, yerr=krat_err, label=name, **kw)

axs2[0].axhline(0, color='gray', lw=0.6, ls='--')
axs2[0].set_xlabel('$y$')
axs2[0].set_ylabel(r'Net-proton $dN/dy$ per event')
axs2[0].set_title(r'Net-Proton ($p_T > 0.4$ GeV/$c$)')
axs2[0].legend(fontsize=7)
axs2[0].axvline(y_beam, color='red', lw=0.7, ls=':', alpha=0.4)
axs2[0].axvline(-y_beam, color='red', lw=0.7, ls=':', alpha=0.4)

axs2[1].axhline(1, color='gray', lw=0.6, ls='--')
axs2[1].set_xlabel('$y$')
axs2[1].set_ylabel(r'$K^+/K^-$')
axs2[1].set_title(r'$K^+/K^-$ ($p_T > 0.2$ GeV/$c$)')
axs2[1].set_ylim(0, 4)
axs2[1].legend(fontsize=7)

fig2.suptitle(r'Au+Au @ 7.7 GeV — Baryon Stopping & Strangeness (with $p_T$ cuts)',
              fontsize=15, fontweight='bold')
fig2.tight_layout()
fig2.savefig(os.path.join(OUTDIR, 'fig2_baryon_stopping.png'), dpi=200, bbox_inches='tight')
print("  Saved fig2_baryon_stopping.png")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 3: dN/dy Ratio — MID-RAPIDITY ONLY, with pT cut
#   |y| < 1.5 for pions/kaons, |y| < 1.0 for protons
#   pT > 0.2 GeV for mesons, pT > 0.4 GeV for baryons
# ═══════════════════════════════════════════════════════════════════════
print("Figure 3: dN/dy ratio (mid-rapidity, pT-cut)...")

# Use finer bins in the restricted range for smooth ratio
ybins_mid = np.arange(-1.5, 1.51, 0.2)  # |y| < 1.5, Δy = 0.2
bc_mid = 0.5 * (ybins_mid[:-1] + ybins_mid[1:])
bw_mid = ybins_mid[1] - ybins_mid[0]

pt_cuts = {
    r'Pions ($\pi^{\pm}$)': 0.2,
    r'Kaons ($K^{\pm}$)': 0.2,
    r'Protons ($p + \bar{p}$)': 0.4,
}

fig3, axs3 = plt.subplots(1, 3, figsize=(18, 5.5))
ref_name = 'Default (No Medium)'

for idx, (sp_label, sp_pids) in enumerate(species_list):
    ax = axs3[idx]
    pt_cut = pt_cuts[sp_label]

    D_ref = DATA[ref_name]
    mask_ref = np.isin(D_ref['pid'], sp_pids) & (D_ref['pt'] > pt_cut)
    h_ref, _ = np.histogram(D_ref['y'][mask_ref], bins=ybins_mid)
    ref_norm = h_ref / D_ref['nev']

    for name in CONFIG_NAMES[1:]:
        if name not in DATA:
            continue
        D = DATA[name]
        mask = np.isin(D['pid'], sp_pids) & (D['pt'] > pt_cut)
        h_mod, _ = np.histogram(D['y'][mask], bins=ybins_mid)
        mod_norm = h_mod / D['nev']

        # Only plot bins with decent stats
        good = (h_ref > 100) & (h_mod > 100)
        ratio = np.full_like(bc_mid, np.nan)
        ratio_err = np.full_like(bc_mid, np.nan)
        ratio[good] = mod_norm[good] / ref_norm[good]
        ratio_err[good] = ratio[good] * np.sqrt(1.0/h_mod[good] + 1.0/h_ref[good])

        ax.errorbar(bc_mid[good], ratio[good], yerr=ratio_err[good],
                    color=COLORS[name], marker=MARKERS[name], ls=LS[name],
                    markersize=5, label=name, capsize=2, elinewidth=0.7)

    ax.axhline(1.0, color='black', lw=1.2, ls='--')
    ax.set_xlabel('Rapidity $y$')
    ax.set_ylabel('Ratio to Default')
    ax.set_title(sp_label + f' ($p_T > {pt_cut}$ GeV/$c$)')
    ax.set_xlim(-1.6, 1.6)
    ax.legend(fontsize=7)

fig3.suptitle(r'Medium Modification Ratio — Mid-Rapidity ($|y|<1.5$), Au+Au @ 7.7 GeV',
              fontsize=15, fontweight='bold')
fig3.tight_layout()
fig3.savefig(os.path.join(OUTDIR, 'dndy_ratio_to_default.png'), dpi=200, bbox_inches='tight')
print("  Saved dndy_ratio_to_default.png")


print("\nDone!")
