import numpy as np

def build_data():
    in_file = 'DataAnalysis/data_fields_18_oct/file_eta0_T0_fs0.txt'
    out_file = 'model_data.csv'
    
    with open(in_file, 'r') as f:
        lines = f.readlines()
        
    out_lines = []
    out_lines.append("density,m_u,m_d,m_s,M_B,V_u,V_d,V_s\n")
    
    for line in lines:
        parts = [p for p in line.strip().split() if p]
        if not parts: continue
        
        rho = float(parts[4])
        m_u = float(parts[15])
        m_d = float(parts[16])
        m_s = float(parts[17])
        
        omega_0 = float(parts[13])
        rho_0 = float(parts[14])
        
        # SU(3) coupling constants from Parameteric.py
        # gv = 10.92
        # gwu = gwd = gru = gv / (2*sqrt(2)) = 3.8608
        g_v = 3.8608 
        
        # Proper CQMF vector potentials from primary meson fields
        v_u = g_v * omega_0 + g_v * rho_0
        v_d = g_v * omega_0 - g_v * rho_0
        v_s = 0.0  # s quark does not couple to omega or rho
        
        # M_B is not really used in ZPC but we can keep a placeholder
        m_b = 939.0
        
        out_lines.append(f"{rho:.4f},{m_u:.4f},{m_d:.4f},{m_s:.4f},{m_b:.4f},{v_u:.4f},{v_d:.4f},{v_s:.4f}\n")
        
    with open(out_file, 'w') as f:
        f.writelines(out_lines)
        
if __name__ == '__main__':
    build_data()
    print("model_data.csv created")
