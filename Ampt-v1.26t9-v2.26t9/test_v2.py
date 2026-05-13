import uproot, numpy as np
for label, fpath in [("Linear", "/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/local_density_approach/ana/ampt_localdensity.root"),
                     ("Gaussian", "/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9/grid_size_10x3/kekcc_submission/ampt_localdensity.root")]:
    d = uproot.open(fpath)["ampt"].arrays(["pid", "px", "py", "pz", "mass"], library="np", entry_stop=1000000)
    e = np.sqrt(d["px"]**2 + d["py"]**2 + d["pz"]**2 + d["mass"]**2)
    y = 0.5 * np.log((e + d["pz"]) / (e - d["pz"]))
    pt = np.sqrt(d["px"]**2 + d["py"]**2)
    phi = np.arctan2(d["py"], d["px"])
    v2 = np.cos(2*phi)
    mask = (np.abs(d["pid"]) == 211) & (np.abs(y) < 1.0) & (pt > 0.5) & (pt < 1.5)
    print(f"{label} Mean v2 (Pions, 0.5<pT<1.5): {np.mean(v2[mask]):.4f}")
