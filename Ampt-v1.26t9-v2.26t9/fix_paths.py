import glob, re
scripts = glob.glob('plot_*.py')
new_dict = """files = {
    'Default': "ana/ampt_default.dat",
    'Fixed Density (1$\\\\rho_0$)': "ana/ampt_modified.dat",
    'Local Density (Phase 1)': "ana/ampt_localdensity.dat"
}"""
for s in scripts:
    with open(s, 'r') as f: content = f.read()
    content = re.sub(r'files\s*=\s*\{[^}]*\}', new_dict, content)
    with open(s, 'w') as f: f.write(content)
print("Updated all python scripts!")
