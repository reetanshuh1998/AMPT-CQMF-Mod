import re

with open('/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/kekcc_20x20x20_production/zpc.f', 'r') as f:
    code = f.read()

# Replace all (10,10,10) arrays in common blocks and declarations
code = code.replace('(10,10,10)', '(20,20,20)')

# Replace loop bounds from 10 to 20
code = re.sub(r'do ix=1,10', 'do ix=1,20', code)
code = re.sub(r'do iy=1,10', 'do iy=1,20', code)
code = re.sub(r'do iz=1,10', 'do iz=1,20', code)

code = re.sub(r'if\(ix\.gt\.10\) ix=10', 'if(ix.gt.20) ix=20', code)
code = re.sub(r'if\(iy\.gt\.10\) iy=10', 'if(iy.gt.20) iy=20', code)
code = re.sub(r'if\(iz\.gt\.10\) iz=10', 'if(iz.gt.20) iz=20', code)

code = re.sub(r'jx\.le\.10', 'jx.le.20', code)
code = re.sub(r'jy\.le\.10', 'jy.le.20', code)
code = re.sub(r'jz\.le\.10', 'jz.le.20', code)

# Replace cell dimension divisions
code = re.sub(r'/10d0', '/20d0', code)

# Add the SAVE statement in build_rhob_grid
save_target = """      integer ix,iy,iz,inside
      save"""
save_replacement = """      integer ix,iy,iz,inside
      save Bcell, Ecell, Pxcell, Pycell, Pzcell
      save"""
code = code.replace(save_target, save_replacement)

with open('/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/kekcc_20x20x20_production/zpc.f', 'w') as f:
    f.write(code)

print("Done fixing zpc.f")
