import numpy as np

def build_data():
    in_file = 'DataAnalysis/data_fields_18_oct/file_eta0_T0_fs0.txt'
    out_file = 'model_data_v2.csv'
    
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
        
        e_u = float(parts[18])
        e_d = float(parts[19])
        e_s = float(parts[20])
        
        v_u = e_u - m_u
        v_d = e_d - m_d
        v_s = e_s - m_s
        
        # M_B is not really used in ZPC but we can keep a placeholder
        m_b = 939.0
        
        out_lines.append(f"{rho:.4f},{m_u:.4f},{m_d:.4f},{m_s:.4f},{m_b:.4f},{v_u:.4f},{v_d:.4f},{v_s:.4f}\n")
        
    with open(out_file, 'w') as f:
        f.writelines(out_lines)
        
if __name__ == '__main__':
    build_data()
    print("model_data_v2.csv created")
