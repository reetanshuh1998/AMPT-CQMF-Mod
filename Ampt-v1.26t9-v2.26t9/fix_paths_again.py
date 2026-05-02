import glob, re
scripts = glob.glob('plot_*.py')
for s in scripts:
    with open(s, 'r') as f: content = f.read()
    content = content.replace("Fixed Density (1$\\rho_0$)", r"Fixed Density ($1\rho_0$)")
    content = content.replace("Fixed Density (1$ho_0$)", r"Fixed Density ($1\rho_0$)")
    with open(s, 'w') as f: f.write(content)
print("Updated all python scripts again!")
