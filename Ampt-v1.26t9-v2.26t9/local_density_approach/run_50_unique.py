import os, subprocess, time, shutil

TARGET_EVENTS = 50
MASTER_FILE = '../ana/ampt_localdensity_fixed.dat'
WORKDIR = 'robust_ld_run'

os.makedirs(WORKDIR, exist_ok=True)
os.chdir(WORKDIR)

for f in ['../ampt', '../model_data.csv', '../input.ampt']:
    if os.path.exists(f):
        shutil.copy(f, '.')

with open('input.ampt', 'r') as f:
    lines = f.readlines()
with open('input.ampt', 'w') as f:
    for line in lines:
        if 'NEVNT' in line:
            f.write('1    ! NEVNT (total number of events)\n')
        elif 'ihjsed' in line and line.strip().startswith('0'):
            f.write('11		! ihjsed: take HIJING seed from below (D=0)or at runtime(11)\n')
        else:
            f.write(line)

with open('input.density', 'w') as f:
    f.write('2.0\n2\n')

os.makedirs('ana', exist_ok=True)

if os.path.exists(MASTER_FILE):
    os.remove(MASTER_FILE)

events_collected = 0
attempts = 0

with open(MASTER_FILE, 'a') as master:
    while events_collected < TARGET_EVENTS:
        attempts += 1
        seed = int(time.time() * 1000) % 1000000 + attempts
        with open('nseed_runtime', 'w') as f:
            f.write(f'{seed}\n')
            
        if os.path.exists('ana/ampt.dat'):
            os.remove('ana/ampt.dat')
            
        try:
            subprocess.run(['./ampt'], stdin=open('nseed_runtime', 'r'), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=180)
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
                    print(f"[{events_collected}/{TARGET_EVENTS}] Event collected.", flush=True)

print("Done!")
