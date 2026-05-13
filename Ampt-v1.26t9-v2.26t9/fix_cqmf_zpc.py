import re

with open('/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/kekcc_20x20x20_production/zpc.f', 'r') as f:
    lines = f.readlines()

out_lines = []
in_cqmf = False

for line in lines:
    if "LOCAL DENSITY CQMF SUBROUTINES" in line:
        in_cqmf = True
    
    if in_cqmf:
        # Replace 10,10,10 with 20,20,20
        line = line.replace('(10,10,10)', '(20,20,20)')
        
        # Replace loop limits
        line = re.sub(r'do (ix|iy|iz)=1,10', r'do \1=1,20', line)
        line = re.sub(r'if\((ix|iy|iz)\.gt\.10\)', r'if(\1.gt.20)', line)
        line = re.sub(r'(ix|iy|iz)=10', r'\1=20', line)
        
        # In smooth_rhob
        line = re.sub(r'(jx|jy|jz)\.le\.10', r'\1.le.20', line)
        
        # In compute_gradients
        line = re.sub(r'(ixp|iyp|izp)=min\((.+?),10\)', r'\1=min(\2,20)', line)
        
        # Cell divisions
        line = re.sub(r'/10d0', '/20d0', line)
        
        # Add SAVE statement in build_rhob_grid
        if "      integer ix,iy,iz,inside" in line:
            out_lines.append(line)
            out_lines.append("      save Bcell, Ecell, Pxcell, Pycell, Pzcell\n")
            continue
            
    out_lines.append(line)

with open('/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/kekcc_20x20x20_production/zpc.f', 'w') as f:
    f.writelines(out_lines)

print("CQMF isolated modifications done.")
