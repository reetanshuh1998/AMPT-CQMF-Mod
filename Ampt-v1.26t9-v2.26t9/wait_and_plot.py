import time, os, subprocess
print("Waiting for ampt_localdensity.dat to be successfully generated...")
while not os.path.exists("local_density_approach/ana/ampt_localdensity.dat"):
    time.sleep(10)
time.sleep(5)  # Wait for file to finish writing
print("Data found! Running plotting scripts...")
subprocess.run("python3 plot_advanced_splittings.py && python3 plot_advanced_baryon_stopping.py && python3 plot_v1_v2.py && python3 plot_pt_spectra.py && python3 plot_proton_kaon_production.py", shell=True)
print("All plots finished successfully!")
