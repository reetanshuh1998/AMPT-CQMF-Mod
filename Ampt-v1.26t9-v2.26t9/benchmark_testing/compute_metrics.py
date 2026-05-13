import numpy as np
import json
import os
import csv

BASE_DIR = '/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9'
DATA_FILE = os.path.join(BASE_DIR, 'benchmark_testing', 'model_v2_data.json')
STAR_DIR = os.path.join(BASE_DIR, 'HEPData-ins1395151-v2-csv')
STAR_SPLIT_P = os.path.join(BASE_DIR, 'star_data', 'v2_splitting_p_pbar_7.7_10_40.csv')
OUT_RESULTS = os.path.join(BASE_DIR, 'benchmark_testing', 'benchmark_results.txt')
OUT_CSV = os.path.join(BASE_DIR, 'benchmark_testing', 'benchmark_results.csv')

def load_star_csv(filename):
    filepath = os.path.join(STAR_DIR, filename)
    pt_list, v2_list, err_list = [], [], []
    with open(filepath) as f:
        for row in csv.reader(f):
            if not row or row[0].startswith('#') or row[0].startswith('PT') or row[0].startswith('$'): continue
            try:
                pt_val, v2_val = float(row[0]), row[1]
                if v2_val == '-': continue
                pt_list.append(pt_val)
                v2_list.append(float(v2_val))
                err_list.append(float(row[2])) # stat error
            except: continue
    return np.array(pt_list), np.array(v2_list), np.array(err_list)

def load_star_splitting(filepath):
    data = np.genfromtxt(filepath, delimiter=',', skip_header=1)
    return data[:, 0], data[:, 1], data[:, 4] # pt, dv2, stat error

star_data = {
    'Pions': load_star_csv('Table110.csv'),
    'Kaons': load_star_csv('Table111.csv'),
    'Protons': load_star_csv('Table107.csv'),
    'Delta_v2': load_star_splitting(STAR_SPLIT_P)
}

def compute_metrics(x_exp, y_exp, err_exp, x_mod, y_mod):
    # Interpolate model to experimental points
    # We only interpolate where model data is valid (not NaN)
    valid_mod = ~np.isnan(y_mod)
    if not np.any(valid_mod):
        return np.nan, np.nan
        
    x_mod_v = np.array(x_mod)[valid_mod]
    y_mod_v = np.array(y_mod)[valid_mod]
    
    y_mod_interp = np.interp(x_exp, x_mod_v, y_mod_v, left=np.nan, right=np.nan)
    
    # Only compute metrics where interpolation is valid
    valid = ~np.isnan(y_mod_interp) & (err_exp > 0)
    
    if np.sum(valid) == 0:
        return np.nan, np.nan
        
    chi2 = np.sum(((y_exp[valid] - y_mod_interp[valid]) / err_exp[valid])**2)
    chi2_ndf = chi2 / np.sum(valid)
    
    rmse = np.sqrt(np.mean((y_exp[valid] - y_mod_interp[valid])**2))
    
    return chi2_ndf, rmse

with open(DATA_FILE, 'r') as f:
    model_data = json.load(f)

results = {}

for model_name, data in model_data.items():
    results[model_name] = {}
    chi2_tot = 0
    
    for species in ['Pions', 'Kaons', 'Protons']:
        x_exp, y_exp, err_exp = star_data[species]
        x_mod = data[species]['pt']
        y_mod = data[species]['v2']
        
        c2, rmse = compute_metrics(x_exp, y_exp, err_exp, x_mod, y_mod)
        results[model_name][species] = {'chi2': c2, 'rmse': rmse}
        if not np.isnan(c2): chi2_tot += c2
        
    # Delta v2
    x_exp, y_exp, err_exp = star_data['Delta_v2']
    
    # Compute model Delta v2
    x_mod = np.array(data['Protons']['pt'])
    y_mod_p = np.array(data['Protons']['v2'])
    y_mod_pb = np.array(data['Antiprotons']['v2'])
    y_mod_dv2 = y_mod_p - y_mod_pb
    
    c2_dv2, rmse_dv2 = compute_metrics(x_exp, y_exp, err_exp, x_mod, y_mod_dv2)
    results[model_name]['Delta_v2'] = {'chi2': c2_dv2, 'rmse': rmse_dv2}
    if not np.isnan(c2_dv2): chi2_tot += c2_dv2
    
    results[model_name]['Global_Chi2'] = chi2_tot

# Output Results
with open(OUT_RESULTS, 'w') as f:
    f.write("Quantitative Benchmark Results: chi^2/ndf (RMSE)\n")
    f.write("="*80 + "\n")
    f.write(f"{'Model':<15} | {'Pions':<12} | {'Kaons':<12} | {'Protons':<12} | {'Delta v2':<12} | {'Global Chi2':<10}\n")
    f.write("-" * 80 + "\n")
    
    for model in sorted(results.keys()):
        r = results[model]
        pi_str = f"{r['Pions']['chi2']:.1f} ({r['Pions']['rmse']:.3f})"
        k_str  = f"{r['Kaons']['chi2']:.1f} ({r['Kaons']['rmse']:.3f})"
        p_str  = f"{r['Protons']['chi2']:.1f} ({r['Protons']['rmse']:.3f})"
        dv_str = f"{r['Delta_v2']['chi2']:.1f} ({r['Delta_v2']['rmse']:.4f})"
        glob   = f"{r['Global_Chi2']:.1f}"
        
        f.write(f"{model:<15} | {pi_str:<12} | {k_str:<12} | {p_str:<12} | {dv_str:<12} | {glob:<10}\n")

with open(OUT_CSV, 'w') as f:
    f.write("Model,Pions_chi2,Pions_rmse,Kaons_chi2,Kaons_rmse,Protons_chi2,Protons_rmse,Deltav2_chi2,Deltav2_rmse,Global_chi2\n")
    for model in sorted(results.keys()):
        r = results[model]
        f.write(f"{model},{r['Pions']['chi2']:.2f},{r['Pions']['rmse']:.4f},{r['Kaons']['chi2']:.2f},{r['Kaons']['rmse']:.4f},{r['Protons']['chi2']:.2f},{r['Protons']['rmse']:.4f},{r['Delta_v2']['chi2']:.2f},{r['Delta_v2']['rmse']:.4f},{r['Global_Chi2']:.2f}\n")

print(f"Results computed and saved to {OUT_RESULTS} and {OUT_CSV}")
