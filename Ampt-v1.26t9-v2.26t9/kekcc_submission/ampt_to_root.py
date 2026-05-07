import ROOT
from array import array
import sys

if len(sys.argv) != 3:
    print("Usage: python3 ampt_to_root.py <input.dat> <output.root>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

f_out = ROOT.TFile(output_file, "RECREATE")
tree = ROOT.TTree("ampt", "AMPT Events")

b = array('f', [0.])
pid = array('i', [0])
px = array('f', [0.])
py = array('f', [0.])
pz = array('f', [0.])
mass = array('f', [0.])

tree.Branch("b", b, "b/F")
tree.Branch("pid", pid, "pid/I")
tree.Branch("px", px, "px/F")
tree.Branch("py", py, "py/F")
tree.Branch("pz", pz, "pz/F")
tree.Branch("mass", mass, "mass/F")

current_b = 0.0
n_left = 0

with open(input_file, 'r') as f_in:
    for line in f_in:
        cols = line.split()
        if not cols: continue
        
        if n_left == 0:
            try:
                n_left = int(cols[2])
                current_b = float(cols[3])
            except:
                pass
        else:
            n_left -= 1
            try:
                pid_val = int(cols[0])
                b[0] = current_b
                pid[0] = pid_val
                px[0] = float(cols[1])
                py[0] = float(cols[2])
                pz[0] = float(cols[3])
                mass[0] = float(cols[4])
                tree.Fill()
            except:
                pass

f_out.Write()
f_out.Close()
print(f"Successfully converted {input_file} to {output_file}")
