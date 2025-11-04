
#%% Imports
import cedalion
import os
import shutil

#%% Load raw files

dyad_num = '04'   # CHANGE
sub_1 = 'nc007' # CHANGE (as seen in file names)
sub_2 = 'pwa006' # CHANGE (as seen in file names)

rename = True
sub_1_rename = 'BUNC277'  #'PWA006'  # CHANGE (check run log)
sub_2_rename = 'BUMA328' # SLP005  # CHANGE (check run log)

# Raw data folder
root_dir = '/projectnb/nphfnirs/s/datasets/Hyperscanning_patients_2025/'
raw_dir = os.path.join(root_dir, 'sourcedata', 'raw', f'dyad_{dyad_num}')

# Open readme file
readme_file = os.path.join(raw_dir, 'README.md')
runs_dict = {}
with open(readme_file, 'r') as file:
    lines = file.read()
    for line in lines.splitlines():  # Changed this line
        # Strip whitespace from the line
        line = line.strip()
        
        # Check if the line matches the pattern "XXX - Name"
        if "=" in line:
            # Split the line into key and value
            key, value = line.split("=", 1)
            # Add the key-value pair to the dictionary
            runs_dict[key.strip()] = value.strip()

# Print the resulting dictionary
print(runs_dict)

#%%Create folder structure for BIDS
runs = list(runs_dict.keys())

# Create new dir for each subject if it doesn't exist
if rename == True:
    sub_dir = os.path.join(root_dir, sub_1_rename)
    sub_dir_2 = os.path.join(root_dir, sub_2_rename)
else:
    sub_dir = os.path.join(root_dir, sub_1)
    sub_dir_2 = os.path.join(root_dir, sub_2)

if not os.path.exists(sub_dir):
    os.makedirs(sub_dir)
                      
if not os.path.exists(os.path.join(sub_dir, f'ses-dyadic{dyad_num}')):
    os.makedirs(os.path.join(sub_dir, f'ses-dyadic{dyad_num}'))
    os.makedirs(os.path.join(sub_dir, f'ses-dyadic{dyad_num}', 'nirs'))

if not os.path.exists(sub_dir_2):
    os.makedirs(sub_dir_2)

if not os.path.exists(os.path.join(sub_dir_2, f'ses-dyadic{dyad_num}')):
    os.makedirs(os.path.join(sub_dir_2, f'ses-dyadic{dyad_num}'))
    os.makedirs(os.path.join(sub_dir_2, f'ses-dyadic{dyad_num}', 'nirs'))

sub_1_nirs_dir = os.path.join(sub_dir, f'ses-dyadic{dyad_num}', 'nirs')
sub_2_nirs_dir = os.path.join(sub_dir_2, f'ses-dyadic{dyad_num}', 'nirs')

#%% Rename files in raw folder

# # Walk through the directory tree
# backup_dir =os.path.join(root_dir, 'sourcedata_backup')
# if 'sourcedata_backup' in raw_dir:
#     backup_dir =os.path.join(root_dir, 'sourcedata_backup_backup')
# if rename == True:
#     if not os.path.exists(os.path.join(backup_dir, 'raw', f'dyad_{dyad_num}')):
#         #os.makedirs(os.path.join(backup_dir, 'raw', f'dyad_{dyad_num}'))
#         shutil.copytree(raw_dir, backup_dir) 
#         print(f"Successfully copied {raw_dir} to {backup_dir}")
#     for root, dirs, files in os.walk(raw_dir, topdown=False):  # Use topdown=False to rename from the bottom up
#         # Rename files
#         for file in files:
#             if file.startswith(sub_1) or file.startswith(sub_2):
#                 old_file_path = os.path.join(root, file)
#                 new_file_name = file.replace(sub_1, sub_1_rename).replace(sub_2, sub_2_rename)
#                 new_file_path = os.path.join(root, new_file_name)
#                 os.rename(old_file_path, new_file_path)
#                 print(f"Renamed file: {old_file_path} -> {new_file_path}")
        
#         # Rename directories
#         for dir in dirs:
#             if dir.startswith(sub_1) or dir.startswith(sub_2):
#                 old_dir_path = os.path.join(root, dir)
#                 new_dir_name = dir.replace(sub_1, sub_1_rename).replace(sub_2, sub_2_rename)
#                 new_dir_path = os.path.join(root, new_dir_name)
#                 os.rename(old_dir_path, new_dir_path)
#                 print(f"Renamed directory: {old_dir_path} -> {new_dir_path}")

#%% Copy snirfs to their respective BIDS folders and rename them 
if rename == True:
    sub_id_1 = sub_1_rename
    sub_id_2 = sub_2_rename
else:
    sub_id_1 = sub_1
    sub_id_2 = sub_2    
 

for run_num in runs_dict:
    run_num_fr = None  # Initialize run_num_fr be none
    run_task = runs_dict[run_num]
    print(run_num)
    print(run_task)
    if 'run' in run_task:
        run_num_fr = run_task.split('run')[1].strip()
        run_task = run_task.split('run')[0].strip()
        if '0' not in run_num_fr:
            run_num_fr = '0' + run_num_fr

    #
    # Copy snirfs to sub_1 nirs dir
    #
    src_snirf = os.path.join(raw_dir, f'{sub_1}_{run_num}', f'{sub_1}_{run_num}.snirf') # get file names
    if run_num_fr:
        dest_snirf = os.path.join(sub_1_nirs_dir, f'{sub_id_1}_ses-dyadic{dyad_num}_task-{run_task}_run-{run_num_fr}_nirs.snirf')
    else:
        dest_snirf = os.path.join(sub_1_nirs_dir, f'{sub_id_1}_ses-dyadic{dyad_num}_task-{run_task}_run-01_nirs.snirf')

    if os.path.exists(src_snirf):   # copy snirfs
        shutil.copy(src_snirf, dest_snirf)
        print(f"Copied {src_snirf} to {dest_snirf}")
    else:
        print(f"Source file {src_snirf} does not exist.")
    
    #
    # Copy snirfs to sub_2
    #
    src_snirf = os.path.join(raw_dir, f'{sub_2}_{run_num}', f'{sub_2}_{run_num}.snirf')
    if run_num_fr:
        dest_snirf = os.path.join(sub_2_nirs_dir, f'{sub_id_2}_ses-dyadic{dyad_num}_task-{run_task}_run-{run_num_fr}_nirs.snirf')
    else:
        dest_snirf = os.path.join(sub_2_nirs_dir, f'{sub_id_2}_ses-dyadic{dyad_num}_task-{run_task}_run-01_nirs.snirf')

    if os.path.exists(src_snirf):
        shutil.copy(src_snirf, dest_snirf)
        print(f"Copied {src_snirf} to {dest_snirf}")
    else:
        print(f"Source file {src_snirf} does not exist.")   





# %%
