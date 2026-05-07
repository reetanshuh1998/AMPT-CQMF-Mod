import os
import subprocess
import time
import shutil

TARGET_EVENTS = 20000

# Use absolute paths so it works no matter where you run it from
BASE_DIR = os.path.abspath('/home/reet/medium_modification/Ampt-v1.26t9-v2.26t9')
WORKDIR = os.path.join(BASE_DIR, 'local_density_approach', 'robust_overnight_run')
MASTER_FILE = os.path.join(BASE_DIR, 'local_density_approach', 'ana', 'ampt_localdensity_fixed.dat')

# Create robust run directory
os.makedirs(WORKDIR, exist_ok=True)
os.chdir(WORKDIR)

# Copy necessary binary and inputs
for f in [os.path.join(BASE_DIR, 'ampt'), 
          os.path.join(BASE_DIR, 'model_data.csv'), 
          os.path.join(BASE_DIR, 'input.ampt')]:
    if os.path.exists(f):
        shutil.copy(f, '.')

# Configure for 1 event per run to bypass crashes safely
with open('input.ampt', 'r') as f:
    lines = f.readlines()
with open('input.ampt', 'w') as f:
    for line in lines:
        if 'NEVNT' in line:
            f.write('1    ! NEVNT (total number of events)\n')
        elif 'ihjsed' in line and line.strip().startswith('0'):
            f.write('11\t\t! ihjsed: take HIJING seed from below (D=0)or at runtime(11)\n')
        else:
            f.write(line)

# Configure local density mode (iqmc=2) and inject the 2.0 HIJING workaround
with open('input.density', 'w') as f:
    f.write('2.0\n2\n')

os.makedirs('ana', exist_ok=True)
os.makedirs(os.path.dirname(MASTER_FILE), exist_ok=True)

# Count already existing events if you are resuming a run
events_collected = 0
if os.path.exists(MASTER_FILE):
    with open(MASTER_FILE, 'r') as master:
        content = master.read()
        events_collected = content.count('      1      1 ') # counts event headers

attempts = 0

print(f"Starting overnight robust data collection.")
print(f"Goal: {TARGET_EVENTS} events.")
print(f"Currently collected: {events_collected} events.")
print(f"Saving safely to: {MASTER_FILE}")

with open(MASTER_FILE, 'a') as master:
    while events_collected < TARGET_EVENTS:
        attempts += 1
        
        # Unique seed generation
        seed = int(time.time() * 1000) % 1000000 + attempts
        with open('nseed_runtime', 'w') as f:
            f.write(f'{seed}\n')
            
        # Clean previous run data
        if os.path.exists('ana/ampt.dat'):
            os.remove('ana/ampt.dat')
            
        # Run AMPT (timeout set to 3 minutes to catch infinite loops)
        try:
            subprocess.run(['./ampt'], stdin=open('nseed_runtime', 'r'), 
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, 
                           timeout=180)
        except subprocess.TimeoutExpired:
            subprocess.run(['pkill', '-f', './ampt'])
            continue
            
        # Check if ampt.dat is generated and non-empty
        if os.path.exists('ana/ampt.dat') and os.path.getsize('ana/ampt.dat') > 0:
            with open('ana/ampt.dat', 'r') as current_data:
                content = current_data.read()
                if '      1      1 ' in content[:50]:  
                    master.write(content)
                    master.flush()
                    events_collected += 1
                    print(f"[{events_collected}/{TARGET_EVENTS}] Event safely appended. (Attempt {attempts})", flush=True)

print("Overnight production successfully hit 2000 events!")
