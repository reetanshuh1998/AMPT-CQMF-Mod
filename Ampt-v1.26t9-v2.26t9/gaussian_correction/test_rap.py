import numpy as np

def rapidity(px, py, pz, m):
    e = np.sqrt(px**2 + py**2 + pz**2 + m**2)
    return 0.5 * np.log((e+pz)/(e-pz+1e-9))

ys = []
with open('scratch_test/ana/ampt.dat', 'r') as f:
    for line in f:
        cols = line.split()
        if len(cols) > 5 and cols[0] == '2212':
            y = rapidity(float(cols[1]), float(cols[2]), float(cols[3]), float(cols[4]))
            ys.append(y)
hist, bins = np.histogram(ys, bins=np.linspace(-4, 4, 30))
print("Proton rapidity distribution:")
print(hist)
