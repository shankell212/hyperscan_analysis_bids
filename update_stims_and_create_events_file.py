#
#
#%% imports
import os
import cedalion

#%% Load in snirf file
root_dir = "/projectnb/nphfnirs/s/datasets/Hyperscanning_controls_2425/"

dyad_num = '11'  # CHANGE
subj = "sub-22"  # CHANGE

ses = 'dyadic' +  dyad_num

tasks = ['IR', 'JR', 'JS','PN'] #, 'IR'] #, 'PN']  # CHANGE
runNums = ['02', '04', '05','02'] #, '02'] #, '02']  # CHANGE

subj_dir = os.path.join(root_dir, subj, f'ses-{ses}', 'nirs')

# Loop through tasks and update stims
for task, runNum in zip(tasks, runNums):
    # load snirf file
    print('task')
    snirf_path = os.path.join(subj_dir, f'{subj}_ses-{ses}_task-{task}_run-{runNum}_nirs.snirf')
    print(f'Loading {snirf_path}')

    # load snirf records
    records = cedalion.io.read_snirf(snirf_path) 
    rec = records[0]

    stim = rec.stim

    sorted_stim = stim.sort_values(by='onset', key=abs)

    if task == 'JS' or task == 'JR' or task == 'IR':
        duration = 120   # 120 sec
    elif task == 'PN':
        duration = 300

    # create a new stim object with the correct duration and onset 
    sorted_stim['trial_type'] = sorted_stim['trial_type'].replace({'1': "rest", '2': task})  # add correct stim names

    # add correct duration
    sorted_stim.loc[sorted_stim['trial_type'] == "rest", 'duration'] = 30.0
    sorted_stim.loc[sorted_stim['trial_type'] == task, 'duration'] = duration

    print(sorted_stim)  # Display the sorted stimulus data


    events_path = os.path.join(subj_dir, f'{subj}_ses-{ses}_task-{task}_run-{runNum}_events.tsv')

    sorted_stim.to_csv(events_path, sep="\t", index=False)



# %%
