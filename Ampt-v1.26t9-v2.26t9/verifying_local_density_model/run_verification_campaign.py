import os
import subprocess
import time
import shutil

TARGET_EVENTS = 20000

BASE_DIR = os.path.abspath('/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9')
WORKDIR = os.path.join(BASE_DIR, 'verifying_local_density_model')

configs = [
    {"name": "ampt_default.dat", "iqmc": 0, "target_rho": 0.0},
    {"name": "ampt_localdensity.dat", "iqmc": 2, "target_rho": 0.0}
]

os.chdir(WORKDIR)

# Copy executable and data
for f in ['ampt', 'model_data.csv']:
    src = os.path.join(BASE_DIR, f)
    if os.path.exists(src):
        shutil.copy(src, '.')

shutil.copy('input.ampt.template', 'input.ampt')
os.makedirs('ana', exist_ok=True)

with open('production.log', 'a') as log:
    log.write(f"--- Starting production campaign: {time.ctime()} ---\n")
    log.write(f"Goal: {TARGET_EVENTS} events per configuration (Min Bias BMAX=15)\n")

for cfg in configs:
    output_file = os.path.join(WORKDIR, 'ana', cfg['name'])
    
    # Write input.density
    with open('input.density', 'w') as f:
        f.write(f"{cfg['target_rho']}\n{cfg['iqmc']}\n")
        
    events_collected = 0
    if os.path.exists(output_file):
        with open(output_file, 'r') as f:
            events_collected = f.read().count('      1      1 ')
            
    attempts = 0
    print(f"\n[{cfg['name']}] Starting production. Current events: {events_collected}/{TARGET_EVENTS}")
    
    with open(output_file, 'a') as master:
        while events_collected < TARGET_EVENTS:
            attempts += 1
            
            # Write runtime seed
            seed = int(time.time() * 1000) % 1000000 + attempts
            with open('nseed_runtime', 'w') as f:
                f.write(f'{seed}\n')
                
            if os.path.exists('ana/ampt.dat'):
                os.remove('ana/ampt.dat')
                
            try:
                subprocess.run(['./ampt'], stdin=open('nseed_runtime', 'r'), 
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, 
                               timeout=120)
            except subprocess.TimeoutExpired:
                subprocess.run(['pkill', '-f', './ampt'])
                continue
                
            if os.path.exists('ana/ampt.dat') and os.path.getsize('ana/ampt.dat') > 0:
                with open('ana/ampt.dat', 'r') as current_data:
                    content = current_data.read()
                    if '      1      1 ' in content[:50]:  
                        master.write(content)
                        master.flush()
                        events_collected += 1
                        if events_collected % 100 == 0:
                            print(f"[{cfg['name']}] Progress: {events_collected}/{TARGET_EVENTS}")
                            with open('production.log', 'a') as log:
                                log.write(f"[{cfg['name']}] Progress: {events_collected}/{TARGET_EVENTS} at {time.ctime()}\n")

    print(f"[{cfg['name']}] Completed {TARGET_EVENTS} events.")
    with open('production.log', 'a') as log:
        log.write(f"[{cfg['name']}] Completed at {time.ctime()}\n")

print("\n--- CAMPAIGN FINISHED ---")
