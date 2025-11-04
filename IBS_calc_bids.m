% Functions for hyperscan analysis
clear
close all
clc
%% IBS calc
% calculate ibs across a single task block and its specific frequency

base_dir = "/projectnb/nphfnirs/s/datasets/Hyperscanning_patients_2025/";
deriv_dir = base_dir + "derivatives";

dyad_num = "03"; %CHANGE FOR EACH DYAD
subject = 'BUMA328'; %CHANGE FOR EACH DYAD, CAN BE EITHER SUBJECT FROM THAT DYAD
task_names = {'JS','IR','JR','PN'}%, 'PN'};
runs = ["01", "01", "01", "01"]%, "02"];

specified_dyad_name = 'dyadic' + dyad_num; %CHANGE FOR EACH DYAD
sub_dir = fullfile(base_dir, subject, "ses-"+specified_dyad_name, "nirs");

%Sets path to allow to get coherence values
results_dir = fullfile(deriv_dir, specified_dyad_name);
addpath(results_dir);

%% Load WTC results
WTC_data = load(char(join(['WTC_GroupRegion_',specified_dyad_name,'.mat'], '')));

fprintf('\n WTC results for dyad %s loaded successfully \n', dyad_num)
%Load individual subject data 
%subject_data = load(char(join([specified_subject,'_data_',specified_dyad_name,'.mat'])));

%%

dyad_data = struct();

%Loop through each task
for i = 1:length(task_names)
%for i = 1
    if contains('04', dyad_num) && i == 5
        calc_block2 = 0;
    else 
        calc_block2 = 1;
    end
    curr_task = task_names{i};
    curr_run = runs(i);

    %task_run = curr_task + "_run-" + curr_run
    task_run = curr_task + "_" + curr_run

    fil_nm = subject+"_ses-"+specified_dyad_name+"_task-"+curr_task+"_run-"+runs(i)+"_events.tsv"; % events file name
    stims = readtable(fullfile(sub_dir, fil_nm), 'FileType','text', 'Delimiter','\t'); % load in stim table
    start_time = stims.onset(strcmp(stims.trial_type, curr_task)); % grab start times for current task 
    end_time = stims.onset(strcmp(stims.trial_type, "rest")); % grab start times for current task 
    task_length = stims.duration(strcmp(stims.trial_type, curr_task));
    task_length = task_length(1); % get duration of task
    
    ind_period = find(WTC_data.WTC.period{i,1} >=3.2 & WTC_data.WTC.period{i,1} <= 12.8);  % find correct period
    
    for j = 1:10
        region_name = sprintf('Region%d',j);
        
        time = WTC_data.WTC.time{i,j};

        % Determine indices that correlate to each stim time
        %[~,ind_stim] = ismember(stim_vals,time);

        [~, ind_stim_start1] = min(abs(time - start_time(1)));   % nearest index in `time`
        if calc_block2
            [~, ind_stim_start2] = min(abs(time - start_time(2)));  % block 2 % nearest index in `time`
        end
        %t_start  = time(idx);
%         if task_run == 'PN_02'
%             [~, ind_stim_end1] = min(abs(time - end_time(1)));   % nearest index in `time`
%             [~, ind_stim_end2] = min(abs(time - end_time(2)));   % nearest index in `time`
%         else
        
        [~, ind_stim_end1] = min(abs(time - end_time(2)));   % block1 nearest index in `time
        if calc_block2
            [~, ind_stim_end2] = min(abs(time - end_time(3)));  % block2 % nearest index in `time`
        end
        
        %Sometimes there is no coherence values, so set the IBS to NaN if
        %no coherence values
        if (isempty(WTC_data.WTC.coherence{i,j}))
            IBS= NaN; 
            block1_avg_coh = NaN;
            if calc_block2
                block2_avg_coh = NaN;
            end
        else
            %Create block for coherence values of each task block for specified
            %period
            block1 = WTC_data.WTC.coherence{i,j}(ind_period,ind_stim_start1:ind_stim_end1);
            if calc_block2
                block2=  WTC_data.WTC.coherence{i,j}(ind_period,ind_stim_start2:ind_stim_end2);
            end
            % calculate avg coherence for each block
            block1_avg_coh = mean(block1(:));
            if calc_block2
                block2_avg_coh = mean(block2(:));
            end

            %calculate IBS
            if calc_block2
                IBS = (block1_avg_coh +block2_avg_coh)/2;
            end
        end
        if calc_block2
            dyad_IBS_results.(task_run).(region_name) = IBS;
            dyad_IBS_block2.(task_run).(region_name) = block2_avg_coh;
        end
        dyad_IBS_block1.(task_run).(region_name) = block1_avg_coh;
        %dyad_IBS_block2.(task_run).(region_name) = block2_avg_coh;
    end

end

%% Save results
save_filepath = results_dir + "/IBS_results_" + specified_dyad_name + ".mat";
save(save_filepath, "dyad_IBS_block1", "dyad_IBS_block2", "dyad_IBS_results")

%%
