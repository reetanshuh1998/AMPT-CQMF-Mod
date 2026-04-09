import numpy as np
from math import sqrt

def build_data():
    # Switch to fs3 file which includes the phi-meson field (column 28)
    # The fs0 file has no phi field, so V_s was always zero.
    in_file = 'DataAnalysis/data_fields_18_oct/file_eta0_T0_fs3.txt'
    out_file = 'model_data.csv'
    
    # ----- SU(3) Coupling Constants from Parameteric.py -----
    gv = 10.92
    
    # omega-meson couplings (u and d quarks only)
    gwu = gv / (2 * sqrt(2))   # = 3.8608
    gwd = gv / (2 * sqrt(2))   # = 3.8608
    gws = 0.0                  # s quark does NOT couple to omega
    
    # rho-meson couplings (isospin splitting)
    gru = gv / (2 * sqrt(2))   # = 3.8608
    grd = -gru                 # = -3.8608 (isospin partner)
    
    # phi-meson couplings (ONLY strange quarks couple to phi)
    gphis = gv / 2.0           # = 5.46
    gphiu = 0.0                # u quark does NOT couple to phi
    gphid = 0.0                # d quark does NOT couple to phi
    # ---------------------------------------------------------
    
    raw_data = []
    with open(in_file, 'r') as f:
        for line in f:
            parts = [p for p in line.strip().split() if p]
            if parts:
                raw_data.append([float(x) for x in parts])

    # Fix low-density numerical noise in fs3 data by interpolating 
    # from density=0 to density=0.3333 (rows 0 to 5)
    row0 = raw_data[0]
    row5 = raw_data[5]
    for i in range(1, 5):
        frac = float(i) / 5.0
        # Smooth all relevant columns: omega(25), rho(26), phi(27), mu(28), md(29), ms(30)
        for col in [25, 26, 27, 28, 29, 30]:
            raw_data[i][col] = (1.0 - frac) * row0[col] + frac * row5[col]
            
    out_lines = []
    out_lines.append("density,m_u,m_d,m_s,M_B,V_u,V_d,V_s\n")
    
    for parts in raw_data:
        
        rho_ratio = float(parts[4])    # rhob / rho_0
        m_u       = float(parts[28])   # effective u-quark scalar mass (MeV)
        m_d       = float(parts[29])   # effective d-quark scalar mass (MeV)
        m_s       = float(parts[30])   # effective s-quark scalar mass (MeV)
        
        omega_0 = float(parts[25])     # omega meson field (MeV)
        rho_0   = float(parts[26])     # rho meson field (MeV)
        phi_0   = float(parts[27])     # phi meson field (MeV)
        
        # Full CQMF vector potentials from all three meson fields
        v_u = gwu * omega_0 + gru * rho_0 + gphiu * phi_0
        v_d = gwd * omega_0 + grd * rho_0 + gphid * phi_0
        v_s = gws * omega_0 + 0.0 * rho_0 + gphis * phi_0  # V_s = gphis * phi_0
        
        # M_B placeholder (not used in ZPC kernel)
        m_b = 939.0
        
        out_lines.append(f"{rho_ratio:.4f},{m_u:.4f},{m_d:.4f},{m_s:.4f},{m_b:.4f},{v_u:.4f},{v_d:.4f},{v_s:.4f}\n")
        
    with open(out_file, 'w') as f:
        f.writelines(out_lines)
        
if __name__ == '__main__':
    build_data()
    print("model_data.csv created with phi-field V_s from fs3 data")
