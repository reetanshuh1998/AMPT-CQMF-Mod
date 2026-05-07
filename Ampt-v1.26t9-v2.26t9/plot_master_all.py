"""
Master Plot Script — AMPT-CQMF, Au+Au @ 7.7 GeV
Regenerates ALL 7 publication plots with:
  - Correct absolute data file paths
  - Per-event normalization (dN/dy per event)
  - All 5 configurations including Local Density
  - Vectorized numpy arrays for extremely fast v2 computation
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys, time

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size':   11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'legend.fontsize': 8.5,
    'lines.linewidth': 1.8,
    'lines.markersize': 4,
})

# ── Absolute paths ──────────────────────────────────────────────────────────
ANA = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/local_density_approach/ana'
OUTDIR = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/publication_plots'
os.makedirs(OUTDIR, exist_ok=True)

CONFIGS = {
    'Default (No Medium)':    f'{ANA}/ampt_default.dat',
    'Fixed ρ=1ρ₀':           f'{ANA}/ampt_fixed_rho1.dat',
    'Fixed ρ=2ρ₀':           f'{ANA}/ampt_fixed_rho2.dat',
    'Fixed ρ=3ρ₀':           f'{ANA}/ampt_fixed_rho3.dat',
    'Local Density (iqmc=2)': f'{ANA}/ampt_localdensity.dat',
}
COLORS  = ['royalblue', 'darkorange', 'forestgreen', 'firebrick', 'purple']
MARKERS = ['o', 's', '^', 'D', 'P']
LS      = ['-', '--', '-.', ':', (0,(3,1,1,1))]

PID = dict(pip=211,pim=-211,kp=321,km=-321,pr=2212,pbar=-2212,la=3122,lb=-3122)

# ── Vectorized Parser ────────────────────────────────────────────────────────
def parse_ampt_vectorized(fpath):
    """
    Reads ampt.dat into flat numpy arrays for fast vectorized operations.
    Returns: nevents, pid, px, py, pz, mass
    """
    if not os.path.exists(fpath):
        print(f'  [MISSING] {fpath}'); return 0, np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
    
    t0 = time.time()
    nevents = 0
    pids, pxs, pys, pzs, ms = [], [], [], [], []
    
    with open(fpath) as f:
        n_left = 0
        for line in f:
            c = line.split()
            if not c: continue
            if n_left == 0:
                try: 
                    n_left = int(c[2])
                    nevents += 1
                except (ValueError, IndexError): pass
            else:
                n_left -= 1
                try:
                    pids.append(int(c[0]))
                    pxs.append(float(c[1]))
                    pys.append(float(c[2]))
                    pzs.append(float(c[3]))
                    ms.append(float(c[4]))
                except (ValueError, IndexError): pass
                
    pids = np.array(pids, dtype=np.int32)
    pxs = np.array(pxs, dtype=np.float32)
    pys = np.array(pys, dtype=np.float32)
    pzs = np.array(pzs, dtype=np.float32)
    ms = np.array(ms, dtype=np.float32)
    
    print(f'  {os.path.basename(fpath)}: {nevents} events, {len(pids)} particles ({time.time()-t0:.1f}s)')
    return nevents, pids, pxs, pys, pzs, ms

# ── Load all data ────────────────────────────────────────────────────────────
print('Parsing all data files (vectorized)...')
import uproot
def parse_ampt_uproot(filepath):
    if not os.path.exists(filepath): return 0, np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
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
        print(f"  {os.path.basename(filepath)} via UPROOT in {time.time()-t0:.2f}s")
        return nevents, pid, px, py, pz, mass
    except Exception as e:
        print(f"Uproot failed: {e}")
        return 0, np.array([]), np.array([]), np.array([]), np.array([]), np.array([])

DATA = {}
for name, fpath in CONFIGS.items():
    print(f'[{name}]')
    if os.path.exists(fpath.replace('.dat', '.root')):
        fpath = fpath.replace('.dat', '.root')
        nev, pids, px, py, pz, m = parse_ampt_uproot(fpath)
    else:
        nev, pids, px, py, pz, m = parse_ampt_vectorized(fpath)
    
    # Pre-compute y, pT, v2 for ALL particles vectorized
    e = np.sqrt(px*px + py*py + pz*pz + m*m)
    denom = e - pz
    valid_y = (denom > 1e-9) & (e > 1e-9)
    y = np.full_like(e, np.nan)
    y[valid_y] = 0.5 * np.log((e[valid_y] + pz[valid_y]) / denom[valid_y])
    
    pT = np.sqrt(px*px + py*py)
    phi = np.arctan2(py, px)
    v2 = np.cos(2*phi)
    
    DATA[name] = {
        'nev': max(nev,1), 
        'pid': pids, 
        'y': y, 
        'pT': pT, 
        'v2': v2
    }

names  = list(CONFIGS.keys())
ybins  = np.linspace(-4, 4, 25);  ybc = .5*(ybins[:-1]+ybins[1:]); ybw = ybins[1]-ybins[0]
ptbins = np.linspace(0.05, 3.0, 30); ptc = .5*(ptbins[:-1]+ptbins[1:]); ptw = ptbins[1]-ptbins[0]
pt_v2bins = np.linspace(0.1, 2.0, 22); ptv2c = .5*(pt_v2bins[:-1]+pt_v2bins[1:])

# ════════════════════════════════════════════════════════════════════════════
# PLOT 1 — pt_spectra_highstats.png   (Kaon & Pion pT, per event)
# ════════════════════════════════════════════════════════════════════════════
print('\nPlot 1: pT spectra...')
fig, (ax1,ax2) = plt.subplots(1,2, figsize=(14,5))
for i,name in enumerate(names):
    D = DATA[name]
    nev, pid, y, pT = D['nev'], D['pid'], D['y'], D['pT']
    kw = dict(color=COLORS[i], marker=MARKERS[i], ls=LS[i], label=name)
    
    for pids_target, ax in [([PID['kp'],PID['km']],ax1), ([PID['pip'],PID['pim']],ax2)]:
        mask = np.isin(pid, pids_target) & (np.abs(y) < 0.5)
        pts = pT[mask]
        if len(pts)==0: continue
        h,_ = np.histogram(pts, bins=ptbins)
        iy  = h / (2*np.pi * ptc * ptw * nev)
        ax.semilogy(ptc, np.where(iy>0,iy,np.nan), **kw)

for ax, title in zip([ax1,ax2],[r'Kaons ($K^{\pm}$, $|y|<0.5$)',
                                r'Pions ($\pi^{\pm}$, $|y|<0.5$)']):
    ax.set_xlabel(r'$p_T$ (GeV/c)'); ax.set_ylabel(r'$(1/2\pi p_T)\,dN/dydp_T$ per event')
    ax.set_title(title); ax.legend()
fig.suptitle(r'$p_T$ Spectra — Au+Au @ 7.7 GeV (normalized per event)', fontsize=14, fontweight='bold')
fig.tight_layout()
fig.savefig(f'{OUTDIR}/pt_spectra_highstats.png', dpi=200, bbox_inches='tight')

# ════════════════════════════════════════════════════════════════════════════
# PLOT 2 — particle_ratios_comparison.png  (K+/K-, p̄/p vs rapidity with pT cuts)
# ════════════════════════════════════════════════════════════════════════════
print('Plot 2: particle ratios...')
fig, (ax1,ax2) = plt.subplots(1,2, figsize=(14,5))
ybins2 = np.linspace(-3,3,20); bc2 = .5*(ybins2[:-1]+ybins2[1:]); bw2 = ybins2[1]-ybins2[0]
for i,name in enumerate(names):
    D = DATA[name]
    pid, y, pT = D['pid'], D['y'], D['pT']
    kw = dict(color=COLORS[i], marker=MARKERS[i], ls=LS[i], label=name)
    
    hkp,_ = np.histogram(y[(pid==PID['kp']) & (pT > 0.2)], bins=ybins2)
    hkm,_ = np.histogram(y[(pid==PID['km']) & (pT > 0.2)], bins=ybins2)
    hpp,_ = np.histogram(y[(pid==PID['pr']) & (pT > 0.4)], bins=ybins2)
    hpb,_ = np.histogram(y[(pid==PID['pbar']) & (pT > 0.4)], bins=ybins2)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        krat = np.where(hkm>0, hkp/hkm, np.nan)
        prat = np.where(hpp>0, hpb/hpp, np.nan)
    ax1.plot(bc2, krat, **kw); ax2.plot(bc2, prat, **kw)

ax1.axhline(1,color='gray',ls='--',lw=0.8); ax1.set_ylim(0,10)
ax1.set_xlabel('Rapidity $y$'); ax1.set_ylabel('$K^+/K^-$'); ax1.set_title('$K^+/K^-$ Ratio ($p_T > 0.2$ GeV/c)'); ax1.legend()
ax2.axhline(1,color='gray',ls='--',lw=0.8); ax2.set_ylim(0,0.02)
ax2.set_xlabel('Rapidity $y$'); ax2.set_ylabel(r'$\bar{p}/p$'); ax2.set_title(r'$\bar{p}/p$ Ratio ($p_T > 0.4$ GeV/c)'); ax2.legend()
fig.suptitle('Particle Ratios — Au+Au @ 7.7 GeV (with $p_T$ cuts)', fontsize=14, fontweight='bold')
fig.tight_layout()
fig.savefig(f'{OUTDIR}/particle_ratios_comparison.png', dpi=200, bbox_inches='tight')

# ════════════════════════════════════════════════════════════════════════════
# PLOT 3 — flow_v1_v2_comparison.png  (v2 for π, K, p vs pT with statistical errors & STAR BES-I Data)
# ════════════════════════════════════════════════════════════════════════════
print('Plot 3: v2 flow...')
fig, axs = plt.subplots(1,3, figsize=(18,5.5))
v2_specs = [([PID['pip'],PID['pim']],'Pions',axs[0]),
            ([PID['kp'],PID['km']], 'Kaons',axs[1]),
            ([PID['pr'],PID['pbar']],'Protons',axs[2])]

def load_star_v2(filename):
    import csv
    filepath = f'/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/HEPData-ins1395151-v2-csv/{filename}'
    pt_list, v2_list, err_list = [], [], []
    with open(filepath, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or row[0].startswith('#') or row[0].startswith('PT') or row[0].startswith('$'):
                continue
            try:
                pt_val = float(row[0])
                v2_val = row[1]
                if v2_val == '-':
                    continue
                v2_val = float(v2_val)
                err_val = float(row[2])
                pt_list.append(pt_val)
                v2_list.append(v2_val)
                err_list.append(err_val)
            except Exception as e:
                continue
    return np.array(pt_list), np.array(v2_list), np.array(err_list)

            
for i,name in enumerate(names):
    D = DATA[name]
    pid, y, pT, v2 = D['pid'], D['y'], D['pT'], D['v2']
    
    for pids_target, label, ax in v2_specs:
        mask = np.isin(pid, pids_target) & (np.abs(y) < 1.0)
        pT_sub = pT[mask]
        v2_sub = v2[mask]
        
        v2_mean = []
        v2_err = []
        for plo, phi in zip(pt_v2bins[:-1], pt_v2bins[1:]):
            in_bin = (pT_sub >= plo) & (pT_sub < phi)
            n_part = np.sum(in_bin)
            if n_part > 5:
                mean_val = np.mean(v2_sub[in_bin])
                std_val = np.std(v2_sub[in_bin])
                v2_mean.append(mean_val)
                v2_err.append(std_val / np.sqrt(n_part))
            else:
                v2_mean.append(np.nan)
                v2_err.append(np.nan)
                
        ax.errorbar(ptv2c, v2_mean, yerr=v2_err,
                    color=COLORS[i], marker=MARKERS[i], ls='none',
                    markersize=5, capsize=1.5, elinewidth=0.6, label=name)

# Overlay STAR BES-I Data (10-40% Centrality)
star_files = {
    'Pions': ('Table110.csv', 'STAR $\pi^+$ BES-I'),
    'Kaons': ('Table111.csv', 'STAR $K^-$ BES-I'),
    'Protons': ('Table107.csv', 'STAR $p$ BES-I')
}
for label, (fn, star_lbl) in star_files.items():
    for _, l, ax in v2_specs:
        if l == label:
            try:
                s_pt, s_v2, s_err = load_star_v2(fn)
                ax.errorbar(s_pt, s_v2, yerr=s_err, fmt='o',
                            mfc='none', mec='black', ecolor='black',
                            markersize=6, capsize=2, elinewidth=0.8,
                            label=star_lbl, zorder=5)
            except Exception as e:
                print(f"Could not load STAR data for {label}: {e}")

for pids_target,label,ax in v2_specs:
    ax.axhline(0,color='gray',ls='--',lw=0.8)
    ax.set_xlabel(r'$p_T$ (GeV/c)'); ax.set_ylabel(r'$v_2 = \langle\cos 2\phi\rangle$')
    ax.set_title(f'$v_2(p_T)$ — {label} (10-40% Centrality)'); ax.legend()
    ax.set_xlim(0.1, 2.0)
fig.suptitle(r'Elliptic Flow $v_2$ vs $p_T$ — AMPT vs STAR BES-I @ 7.7 GeV', fontsize=14, fontweight='bold')
fig.tight_layout()
fig.savefig(f'{OUTDIR}/flow_v1_v2_comparison.png', dpi=200, bbox_inches='tight')

# ════════════════════════════════════════════════════════════════════════════
# PLOT 4 — baryon_stopping_comparison.png (net-p dN/dy, K+/K-)
# ════════════════════════════════════════════════════════════════════════════
print('Plot 4: baryon stopping...')
fig, axs = plt.subplots(1,3, figsize=(18,5))
for i,name in enumerate(names):
    D = DATA[name]
    nev, pid, y = D['nev'], D['pid'], D['y']
    kw = dict(color=COLORS[i], marker=MARKERS[i], ls=LS[i], label=name)
    
    hp,_  = np.histogram(y[pid==PID['pr']],   bins=ybins)
    hpb,_ = np.histogram(y[pid==PID['pbar']], bins=ybins)
    hkp,_ = np.histogram(y[pid==PID['kp']],   bins=ybins)
    hkm,_ = np.histogram(y[pid==PID['km']],   bins=ybins)
    hla,_ = np.histogram(y[pid==PID['la']],   bins=ybins)
    hlb,_ = np.histogram(y[pid==PID['lb']],   bins=ybins)
    
    netp = (hp - hpb) / (nev*ybw)
    axs[0].plot(ybc, netp, **kw)
    with np.errstate(divide='ignore', invalid='ignore'):
        krat = np.where(hkm>0, hkp/hkm, np.nan)
        lrat = np.where(hlb>0, hla/hlb, np.nan)
    axs[1].plot(ybc, krat, **kw)
    axs[2].plot(ybc, lrat, **kw)

axs[0].axhline(0,color='gray',ls='--',lw=0.8)
axs[0].set_xlabel('$y$'); axs[0].set_ylabel(r'Net-$p$ $dN/dy$ per event')
axs[0].set_title('Net-Proton Baryon Stopping'); axs[0].legend()
for ax,title,ylim in zip(axs[1:],['$K^+/K^-$ vs $y$',r'$\Lambda/\bar{\Lambda}$ vs $y$'],[(0,5),(0,10)]):
    ax.axhline(1,color='gray',ls='--',lw=0.8)
    ax.set_xlabel('$y$'); ax.set_ylim(*ylim); ax.set_title(title); ax.legend()
fig.suptitle('Baryon Stopping & Strangeness — Au+Au @ 7.7 GeV', fontsize=14, fontweight='bold')
fig.tight_layout()
fig.savefig(f'{OUTDIR}/baryon_stopping_comparison.png', dpi=200, bbox_inches='tight')

# ════════════════════════════════════════════════════════════════════════════
# PLOT 5 — proton_v2_splitting_comparison.png (p vs p̄ v2 with statistical errors)
# ════════════════════════════════════════════════════════════════════════════
print('Plot 5: p vs pbar v2 splitting...')
fig, axs = plt.subplots(1,2, figsize=(13,5.5))
for i,name in enumerate(names):
    D = DATA[name]
    pid, y, pT, v2 = D['pid'], D['y'], D['pT'], D['v2']
    
    v2p_val, v2pb_val = [], []
    errp_val, errpb_val = [], []
    
    for pids_target, mkr, lbl, v2_out, err_out in [([PID['pr']], 'o', 'p', v2p_val, errp_val),
                                                   ([PID['pbar']], 's', 'p̄', v2pb_val, errpb_val)]:
        mask = np.isin(pid, pids_target) & (np.abs(y) < 1.0)
        pT_sub = pT[mask]
        v2_sub = v2[mask]
        
        v2_mean = []
        v2_err = []
        for plo, phi in zip(pt_v2bins[:-1], pt_v2bins[1:]):
            in_bin = (pT_sub >= plo) & (pT_sub < phi)
            n_part = np.sum(in_bin)
            if n_part > 3:
                mean_val = np.mean(v2_sub[in_bin])
                std_val = np.std(v2_sub[in_bin])
                v2_mean.append(mean_val)
                v2_err.append(std_val / np.sqrt(n_part))
            else:
                v2_mean.append(np.nan)
                v2_err.append(np.nan)
        
        axs[0].errorbar(ptv2c, v2_mean, yerr=v2_err,
                        color=COLORS[i], marker=mkr, ls='none',
                        markersize=4.5, capsize=1.5, elinewidth=0.6,
                        label=f'{name} ({lbl})')
        v2_out.extend(v2_mean)
        err_out.extend(v2_err)
        
    split = np.array(v2p_val) - np.array(v2pb_val)
    split_err = np.sqrt(np.array(errp_val)**2 + np.array(errpb_val)**2)
    axs[1].errorbar(ptv2c, split, yerr=split_err,
                    color=COLORS[i], marker=MARKERS[i], ls='none',
                    markersize=5, capsize=1.5, elinewidth=0.6, label=name)

axs[0].axhline(0,color='gray',ls='--',lw=0.8)
axs[0].set_xlabel(r'$p_T$ (GeV/c)'); axs[0].set_ylabel(r'$v_2$')
axs[0].set_title(r'$v_2(p_T)$: Proton vs Antiproton')
axs[0].set_xlim(0.1, 2.0)
axs[0].legend(fontsize=7, loc='upper left')

axs[1].axhline(0,color='gray',ls='--',lw=0.8)
axs[1].set_xlabel(r'$p_T$ (GeV/c)'); axs[1].set_ylabel(r'$\Delta v_2 = v_2(p) - v_2(\bar{p})$')
axs[1].set_title(r'$v_2$ Splitting: Sensitivity to Vector Potential')
axs[1].set_xlim(0.1, 2.0)
axs[1].legend(fontsize=8, loc='upper right')

fig.suptitle(r'Proton/$\bar{p}$ $v_2$ Splitting — Au+Au @ 7.7 GeV (Experimental Standard)', fontsize=14, fontweight='bold')
fig.tight_layout()
fig.savefig(f'{OUTDIR}/proton_v2_splitting_comparison.png', dpi=200, bbox_inches='tight')

# ════════════════════════════════════════════════════════════════════════════
# PLOT 6 — proton_kaon_production_highstats.png (dN/dy for K,p,Λ)
# ════════════════════════════════════════════════════════════════════════════
print('Plot 6: proton/kaon production...')
fig, axs = plt.subplots(1,3, figsize=(18,5))
specs = [([PID['kp'],PID['km']],'$K^{\pm}$',axs[0]),
         ([PID['pr'],PID['pbar']],'$p$, $\\bar{p}$',axs[1]),
         ([PID['la'],PID['lb']],'$\\Lambda$, $\\bar{\\Lambda}$',axs[2])]
for i,name in enumerate(names):
    D = DATA[name]
    nev, pid, y = D['nev'], D['pid'], D['y']
    kw = dict(color=COLORS[i], marker=MARKERS[i], ls=LS[i], label=name)
    for pids_target, label, ax in specs:
        mask = np.isin(pid, pids_target)
        h,_ = np.histogram(y[mask], bins=ybins)
        ax.plot(ybc, h/(nev*ybw), **kw)
for pids_target,label,ax in specs:
    ax.set_xlabel('Rapidity $y$'); ax.set_ylabel('$dN/dy$ per event')
    ax.set_title(label+' Rapidity Dist.'); ax.set_ylim(bottom=0); ax.legend()
fig.suptitle('Particle Production — Au+Au @ 7.7 GeV (per event)', fontsize=14, fontweight='bold')
fig.tight_layout()
fig.savefig(f'{OUTDIR}/proton_kaon_production_highstats.png', dpi=200, bbox_inches='tight')

# ════════════════════════════════════════════════════════════════════════════
# PLOT 7 — advanced_observables_highstats.png (mean pT, total yields, K/π)
# ════════════════════════════════════════════════════════════════════════════
print('Plot 7: advanced observables...')
fig, axs = plt.subplots(1,3, figsize=(18,5))

mean_pt_k, mean_pt_pi, mean_pt_p = [], [], []
kpi_ratio, yield_k, yield_pi     = [], [], []

for name in names:
    D = DATA[name]
    nev, pid, y, pT = D['nev'], D['pid'], D['y'], D['pT']
    
    mask_k  = np.isin(pid, [PID['kp'],PID['km']]) & (np.abs(y)<0.5)
    mask_pi = np.isin(pid, [PID['pip'],PID['pim']]) & (np.abs(y)<0.5)
    mask_p  = np.isin(pid, [PID['pr'],PID['pbar']]) & (np.abs(y)<0.5)
    
    mean_pt_k.append( np.mean(pT[mask_k])  if np.sum(mask_k)>0  else 0)
    mean_pt_pi.append(np.mean(pT[mask_pi]) if np.sum(mask_pi)>0 else 0)
    mean_pt_p.append( np.mean(pT[mask_p])  if np.sum(mask_p)>0  else 0)
    
    yk  = np.sum(np.isin(pid, [PID['kp'],PID['km']])) / nev
    ypi = np.sum(np.isin(pid, [PID['pip'],PID['pim']])) / nev
    yield_k.append(yk); yield_pi.append(ypi)
    kpi_ratio.append(yk/ypi if ypi>0 else 0)

x = np.arange(len(names)); w=0.55
bars = axs[0].bar(x, mean_pt_k,  w, color=COLORS[:len(names)], alpha=0.7)
axs[0].set_xticks(x); axs[0].set_xticklabels(names, rotation=20, ha='right', fontsize=8)
axs[0].set_ylabel(r'$\langle p_T\rangle$ (GeV/c)'); axs[0].set_title(r'Mean $p_T$: Kaons ($|y|<0.5$)')

axs[1].bar(x, yield_k,  w*0.4, color=COLORS[:len(names)], alpha=0.9, label='Kaons')
axs[1].bar(x+w*0.4, yield_pi, w*0.4, color=COLORS[:len(names)], alpha=0.5, label='Pions')
axs[1].set_xticks(x+w*0.2); axs[1].set_xticklabels(names, rotation=20, ha='right', fontsize=8)
axs[1].set_ylabel('$dN/dy$ per event'); axs[1].set_title('Total Yields per Event')
axs[1].legend()

axs[2].bar(x, kpi_ratio, w, color=COLORS[:len(names)], alpha=0.8)
axs[2].set_xticks(x); axs[2].set_xticklabels(names, rotation=20, ha='right', fontsize=8)
axs[2].set_ylabel('$K/\\pi$ ratio'); axs[2].set_title('Kaon-to-Pion Ratio')

fig.suptitle('Global Observables — Au+Au @ 7.7 GeV', fontsize=14, fontweight='bold')
fig.tight_layout()
fig.savefig(f'{OUTDIR}/advanced_observables_highstats.png', dpi=200, bbox_inches='tight')

print('\n✓ All 7 plots regenerated successfully (vectorized).')
print(f'  Output directory: {OUTDIR}')
