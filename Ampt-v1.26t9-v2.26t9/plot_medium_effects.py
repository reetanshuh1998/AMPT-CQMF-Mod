"""
AMPT-CQMF Medium Effect Plots
Au+Au @ 7.7 GeV — Plots that directly demonstrate the medium modification

Generates:
  Plot 1: pT spectrum ratio (Modified / Default) — the "money plot"
  Plot 2: Mean pT comparison across models
  Plot 3: Invariant yield spectra (log scale) — slope comparison
  Plot 4: Inverse slope parameter T_eff from Boltzmann fit
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os, time
import uproot

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 11,
    'axes.labelsize': 13,
    'axes.titlesize': 14,
    'legend.fontsize': 9,
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

# ── Parse all ROOT files ─────────────────────────────────────────────────
def parse_ampt_uproot(filepath):
    if not os.path.exists(filepath):
        return 0, None, None, None, None, None
    t0 = time.time()
    try:
        with uproot.open(filepath) as f:
            tree = f["ampt"]
            data = tree.arrays(["b", "pid", "px", "py", "pz", "mass"], library="np")
        b = data["b"]
        nevents = len(np.unique(b))
        pid = data["pid"]
        px = data["px"]
        py = data["py"]
        pz = data["pz"]
        mass = data["mass"]

        pt = np.sqrt(px**2 + py**2)
        e = np.sqrt(px**2 + py**2 + pz**2 + mass**2)
        denom = e - pz
        valid = (denom > 1e-9) & (e > 1e-9)
        y = np.full_like(e, np.nan)
        y[valid] = 0.5 * np.log((e[valid] + pz[valid]) / denom[valid])

        print(f"  Parsed {os.path.basename(filepath)} via UPROOT: "
              f"{nevents} events, {len(pid)} particles ({time.time()-t0:.1f}s)")
        return nevents, pid, pt, y, px, py
    except Exception as ex:
        print(f"  Uproot failed on {filepath}: {ex}")
        return 0, None, None, None, None, None

print("=" * 60)
print("Parsing all ROOT files...")
print("=" * 60)
DATA = {}
for name, fpath in CONFIGS.items():
    print(f"\n[{name}]")
    nev, pid, pt, y, px, py = parse_ampt_uproot(fpath)
    if nev > 0:
        DATA[name] = {'nev': nev, 'pid': pid, 'pt': pt, 'y': y, 'px': px, 'py': py}
print("\n" + "=" * 60)
print("All data parsed. Generating plots...")
print("=" * 60)

# ── Binning setup ────────────────────────────────────────────────────────
pt_bins = np.linspace(0.05, 3.0, 30)
ptc = 0.5 * (pt_bins[:-1] + pt_bins[1:])
ptw = pt_bins[1] - pt_bins[0]

SPECIES = {
    r'Pions ($\pi^{\pm}$)': [PID['pip'], PID['pim']],
    r'Kaons ($K^{\pm}$)': [PID['kp'], PID['km']],
    r'Protons ($p, \bar{p}$)': [PID['pr'], PID['pbar']],
}

# Pre-compute pT histograms for all models/species (|y| < 0.5)
HISTS = {}  # (name, species_label) -> (hist_counts, nev)
for name, D in DATA.items():
    nev = D['nev']
    pid = D['pid']
    pt = D['pt']
    y = D['y']
    mid_rap = np.abs(y) < 0.5
    for sp_label, sp_pids in SPECIES.items():
        mask = np.isin(pid, sp_pids) & mid_rap
        h, _ = np.histogram(pt[mask], bins=pt_bins)
        HISTS[(name, sp_label)] = (h, nev)


# ════════════════════════════════════════════════════════════════════════
# PLOT 1: pT Spectrum Ratio (Modified / Default) — THE MONEY PLOT
# ════════════════════════════════════════════════════════════════════════
print("\nPlot 1: pT spectrum ratio (Modified / Default)...")

fig1, axs1 = plt.subplots(1, 3, figsize=(18, 5.5))
ref_name = 'Default (No Medium)'

for idx, (sp_label, sp_pids) in enumerate(SPECIES.items()):
    ax = axs1[idx]
    ref_h, ref_nev = HISTS[(ref_name, sp_label)]
    ref_norm = ref_h / ref_nev  # per-event yield

    for name in CONFIG_NAMES[1:]:  # skip Default
        if name not in DATA:
            continue
        h, nev = HISTS[(name, sp_label)]
        mod_norm = h / nev

        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.where(ref_norm > 0, mod_norm / ref_norm, np.nan)

        ax.plot(ptc, ratio, color=COLORS[name], marker=MARKERS[name],
                ls=LS[name], markersize=5, label=name)

    ax.axhline(1.0, color='black', lw=1.2, ls='--')
    ax.set_xlabel(r'$p_T$ (GeV/$c$)')
    ax.set_ylabel('Yield Ratio (Modified / Default)')
    ax.set_title(sp_label + r', $|y| < 0.5$')
    ax.set_ylim(0.5, 1.8)
    ax.set_xlim(0.0, 2.5)
    ax.legend(fontsize=8)
    ax.fill_between([0, 0.5], 0.5, 1.8, alpha=0.05, color='green', label='_nolegend_')
    ax.fill_between([1.5, 2.5], 0.5, 1.8, alpha=0.05, color='red', label='_nolegend_')
    ax.text(0.15, 1.65, 'soft\nenhancement', fontsize=7, color='green', ha='center', style='italic')
    ax.text(2.1, 0.6, 'high-$p_T$\nsuppression', fontsize=7, color='red', ha='center', style='italic')

fig1.suptitle(r'Medium Modification of $p_T$ Spectra — Au+Au @ 7.7 GeV ($|y|<0.5$)',
              fontsize=15, fontweight='bold')
fig1.tight_layout()
out1 = os.path.join(OUTDIR, 'medium_pt_ratio.png')
fig1.savefig(out1, dpi=200, bbox_inches='tight')
print(f"  Saved: {out1}")


# ════════════════════════════════════════════════════════════════════════
# PLOT 2: Mean pT Comparison
# ════════════════════════════════════════════════════════════════════════
print("\nPlot 2: Mean pT comparison...")

fig2, axs2 = plt.subplots(1, 3, figsize=(18, 5.5))

for idx, (sp_label, sp_pids) in enumerate(SPECIES.items()):
    ax = axs2[idx]
    means = []
    errs = []
    labels = []
    colors = []

    for name in CONFIG_NAMES:
        if name not in DATA:
            continue
        D = DATA[name]
        mask = np.isin(D['pid'], sp_pids) & (np.abs(D['y']) < 0.5) & (D['pt'] > 0.1)
        pts = D['pt'][mask]
        if len(pts) > 100:
            means.append(np.mean(pts))
            errs.append(np.std(pts) / np.sqrt(len(pts)))
            labels.append(name)
            colors.append(COLORS[name])

    x = np.arange(len(labels))
    bars = ax.bar(x, means, 0.6, yerr=errs, color=colors, alpha=0.8,
                  edgecolor='black', linewidth=0.5, capsize=4)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha='right', fontsize=7)
    ax.set_ylabel(r'$\langle p_T \rangle$ (GeV/$c$)')
    ax.set_title(sp_label + r', $|y|<0.5$')

    # Add value labels on bars
    for bar, val in zip(bars, means):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f'{val:.3f}', ha='center', va='bottom', fontsize=7)

fig2.suptitle(r'Mean Transverse Momentum $\langle p_T \rangle$ — Au+Au @ 7.7 GeV',
              fontsize=15, fontweight='bold')
fig2.tight_layout()
out2 = os.path.join(OUTDIR, 'mean_pt_comparison.png')
fig2.savefig(out2, dpi=200, bbox_inches='tight')
print(f"  Saved: {out2}")


# ════════════════════════════════════════════════════════════════════════
# PLOT 3: Invariant Yield Spectra (log scale — slope comparison)
# ════════════════════════════════════════════════════════════════════════
print("\nPlot 3: Invariant yield spectra...")

fig3, axs3 = plt.subplots(1, 3, figsize=(18, 5.5))

for idx, (sp_label, sp_pids) in enumerate(SPECIES.items()):
    ax = axs3[idx]
    for name in CONFIG_NAMES:
        if name not in DATA:
            continue
        h, nev = HISTS[(name, sp_label)]
        # Invariant yield: (1/2π pT) dN/(dy dpT) per event
        inv_yield = h / (2 * np.pi * ptc * ptw * nev)
        valid = inv_yield > 0
        ax.semilogy(ptc[valid], inv_yield[valid], color=COLORS[name],
                    marker=MARKERS[name], ls=LS[name], markersize=4, label=name)

    ax.set_xlabel(r'$p_T$ (GeV/$c$)')
    ax.set_ylabel(r'$\frac{1}{2\pi p_T} \frac{d^2N}{dydp_T}$ per event')
    ax.set_title(sp_label + r', $|y|<0.5$')
    ax.legend(fontsize=7)
    ax.set_xlim(0, 2.5)

fig3.suptitle(r'Invariant $p_T$ Spectra — Au+Au @ 7.7 GeV ($|y|<0.5$)',
              fontsize=15, fontweight='bold')
fig3.tight_layout()
out3 = os.path.join(OUTDIR, 'invariant_yield_spectra.png')
fig3.savefig(out3, dpi=200, bbox_inches='tight')
print(f"  Saved: {out3}")


# ════════════════════════════════════════════════════════════════════════
# PLOT 4: Inverse Slope T_eff from Boltzmann Fit
# ════════════════════════════════════════════════════════════════════════
print("\nPlot 4: Inverse slope parameter T_eff...")

def boltzmann_mt(mt, A, T):
    """Boltzmann thermal model: dN/dmT ~ A * mT * exp(-mT/T)"""
    return A * mt * np.exp(-mt / T)

# Masses for mT calculation
MASSES = {
    r'Pions ($\pi^{\pm}$)': 0.1396,
    r'Kaons ($K^{\pm}$)': 0.4937,
    r'Protons ($p, \bar{p}$)': 0.9383,
}

fig4, axs4 = plt.subplots(1, 3, figsize=(18, 5.5))

for idx, (sp_label, sp_pids) in enumerate(SPECIES.items()):
    ax = axs4[idx]
    m0 = MASSES[sp_label]
    temps = []
    temp_errs = []
    labels = []
    colors = []

    for name in CONFIG_NAMES:
        if name not in DATA:
            continue
        h, nev = HISTS[(name, sp_label)]
        # Convert pT bins to mT
        mt = np.sqrt(ptc**2 + m0**2)
        inv_yield = h / (2 * np.pi * ptc * ptw * nev)

        # Fit in range where we have good stats (mT - m0 < 1.5 GeV)
        fit_mask = (inv_yield > 0) & ((mt - m0) < 1.5) & ((mt - m0) > 0.05)
        if np.sum(fit_mask) < 4:
            continue

        try:
            popt, pcov = curve_fit(boltzmann_mt, mt[fit_mask], inv_yield[fit_mask],
                                   p0=[1e3, 0.2], maxfev=5000)
            T_fit = popt[1] * 1000  # Convert to MeV
            T_err = np.sqrt(pcov[1, 1]) * 1000
            temps.append(T_fit)
            temp_errs.append(T_err)
            labels.append(name)
            colors.append(COLORS[name])
        except Exception:
            pass

    if len(temps) > 0:
        x = np.arange(len(labels))
        bars = ax.bar(x, temps, 0.6, yerr=temp_errs, color=colors, alpha=0.8,
                      edgecolor='black', linewidth=0.5, capsize=4)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25, ha='right', fontsize=7)
        ax.set_ylabel(r'$T_{\mathrm{eff}}$ (MeV)')
        ax.set_title(sp_label)

        for bar, val in zip(bars, temps):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                    f'{val:.0f}', ha='center', va='bottom', fontsize=8)

fig4.suptitle(r'Inverse Slope Parameter $T_{\mathrm{eff}}$ from Boltzmann Fit — Au+Au @ 7.7 GeV',
              fontsize=15, fontweight='bold')
fig4.tight_layout()
out4 = os.path.join(OUTDIR, 'inverse_slope_Teff.png')
fig4.savefig(out4, dpi=200, bbox_inches='tight')
print(f"  Saved: {out4}")


print("\n" + "=" * 60)
print("All 4 medium-effect plots generated successfully!")
print(f"Output directory: {OUTDIR}")
print("=" * 60)
