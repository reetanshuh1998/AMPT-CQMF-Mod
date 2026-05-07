import os
import subprocess
import time
import shutil

TARGET_EVENTS = 50
MASTER_FILE = '../ana/ampt_localdensity_fixed.dat'
WORKDIR = 'robust_ld_run'

# Create robust run directory
os.makedirs(WORKDIR, exist_ok=True)
os.chdir(WORKDIR)

# Copy necessary files
for f in ['../ampt', '../model_data.csv', '../input.ampt']:
    if os.path.exists(f):
        shutil.copy(f, '.')

# Configure for 1 event
with open('input.ampt', 'r') as f:
    lines = f.readlines()
with open('input.ampt', 'w') as f:
    for line in lines:
        if 'NEVNT' in line:
            f.write('1    ! NEVNT (total number of events)\n')
        else:
            f.write(line)

# Configure local density mode (iqmc=2)
with open('input.density', 'w') as f:
    f.write('2.0\n2\n')

# Ensure output directory exists
os.makedirs('ana', exist_ok=True)

# Clear old master file if starting fresh
if os.path.exists(MASTER_FILE):
    os.remove(MASTER_FILE)

events_collected = 0
attempts = 0

print(f"Starting robust event collection for Local Density Mode. Target: {TARGET_EVENTS} events.")

with open(MASTER_FILE, 'w') as master:
    while events_collected < TARGET_EVENTS:
        attempts += 1
        # Seed generation
        seed = int(time.time() * 1000) % 1000000 + attempts
        with open('nseed_runtime', 'w') as f:
            f.write(f'{seed}\n')
            
        # Clean previous run data
        if os.path.exists('ana/ampt.dat'):
            os.remove('ana/ampt.dat')
            
        # Run AMPT
        try:
            subprocess.run(['./ampt'], stdin=open('nseed_runtime', 'r'), 
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, 
                           timeout=180) # 60 sec timeout in case of infinite loop
        except subprocess.TimeoutExpired:
            # Infinite loop detected, kill it and retry
            subprocess.run(['pkill', '-f', './ampt'])
            continue
            
        # Check if ampt.dat is generated and non-empty
        if os.path.exists('ana/ampt.dat') and os.path.getsize('ana/ampt.dat') > 0:
            with open('ana/ampt.dat', 'r') as current_data:
                # Read to verify it contains event header
                content = current_data.read()
                if '      1      1 ' in content[:50]:  # ievt=1, irun=1 is typical header
                    # Adjust event number so plots know total events
                    # The header line looks like: "      1      1  1025   2.28093..."
                    # We can just write the content directly. The plotting script parses events by counting headers!
                    master.write(content)
                    master.flush()
                    events_collected += 1
                    print(f"[{events_collected}/{TARGET_EVENTS}] Event collected successfully. (Attempt {attempts})")
        else:
            # Failed to generate data (HIJING crash)
            pass

print("Finished data collection successfully!")
