import csv
import os

def parse_hepdata_csv_robust(filename):
    """
    Parses a HEPData CSV file manually.
    Blocks are separated by lines starting with '#: SQRT(S)' or empty lines.
    """
    blocks = []
    current_data = []
    current_header = None
    
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                if current_data:
                    blocks.append((current_header, current_data))
                    current_data = []
                    current_header = None
                continue
            
            if row[0].startswith('#'):
                # Check if this is a new block metadata
                if 'SQRT(S)' in row[0] and current_data:
                    blocks.append((current_header, current_data))
                    current_data = []
                    current_header = None
                continue
                
            if current_header is None:
                current_header = row
            else:
                # Keep all data, even with '-'
                current_data.append(row)
                
    if current_data:
        blocks.append((current_header, current_data))
        
    return blocks

def write_clean_csv(filename, header, data):
    # Filter out rows with '-' in the v2 column (index 1)
    clean_data = [row for row in data if len(row) > 1 and row[1] != '-']
    if clean_data:
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(clean_data)

os.makedirs('star_data', exist_ok=True)

# 1. Protons and Antiprotons (10-40%) @ 7.7 GeV
p_blocks = parse_hepdata_csv_robust('HEPData-ins1210464-v1-csv/Table216.csv')
if len(p_blocks) >= 2:
    write_clean_csv('star_data/v2_proton_7.7_10_40.csv', p_blocks[0][0], p_blocks[0][1])
    write_clean_csv('star_data/v2_pbar_7.7_10_40.csv', p_blocks[1][0], p_blocks[1][1])

# 1b. Protons and Antiprotons (0-10%) @ 7.7 GeV
p010_blocks = parse_hepdata_csv_robust('HEPData-ins1210464-v1-csv/Table204.csv')
if len(p010_blocks) >= 2:
    write_clean_csv('star_data/v2_proton_7.7_0_10.csv', p010_blocks[0][0], p010_blocks[0][1])
    write_clean_csv('star_data/v2_pbar_7.7_0_10.csv', p010_blocks[1][0], p010_blocks[1][1])

# 2. Charged Kaons (0-80%) @ 7.7 GeV
k_blocks = parse_hepdata_csv_robust('HEPData-ins1210464-v1-csv/Table133.csv')
if len(k_blocks) >= 2:
    write_clean_csv('star_data/v2_kp_7.7_0_80.csv', k_blocks[0][0], k_blocks[0][1])
    write_clean_csv('star_data/v2_km_7.7_0_80.csv', k_blocks[1][0], k_blocks[1][1])

# 3. Pions (0-80%) @ 7.7 GeV
pi_blocks = parse_hepdata_csv_robust('HEPData-ins1210464-v1-csv/Table1.csv')
if len(pi_blocks) >= 1:
    write_clean_csv('star_data/v2_pip_7.7_0_80.csv', pi_blocks[0][0], pi_blocks[0][1])

# 4. Delta v2 (p - pbar) (10-40%) @ 7.7 GeV
dv2_blocks = parse_hepdata_csv_robust('HEPData-ins1210464-v1-csv/Table217.csv')
if len(dv2_blocks) >= 1:
    write_clean_csv('star_data/v2_splitting_p_pbar_7.7_10_40.csv', dv2_blocks[0][0], dv2_blocks[0][1])

print("STAR data extraction complete (robust mode). Files saved in star_data/")
