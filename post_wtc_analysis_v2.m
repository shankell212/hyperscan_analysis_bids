clear
clc
close all

%%
% Analyze hyperscanning data after WTC

base_dir = "/projectnb/nphfnirs/s/datasets/Hyperscanning_patients_2025/";
deriv_dir = base_dir + "derivatives/";

save_dir = deriv_dir + "plots/";

filesAndDirs = dir(deriv_dir);
names = string({filesAndDirs.name});
%disp(names)

dyad_lst = names(startsWith(names, "dyadic"));

task_lst = ["JS_01", "IR_01", "JR_01", "PN_01"];

%% Load in excel run log
opts = detectImportOptions(base_dir + "Hyperscanning_Log.xlsx",'FileType', 'text', 'Delimiter', ',');  % force CSV/text
opts.VariableNamesLine = 2;  % headers on line 2
opts.DataLines         = [3 Inf];  % data starts on row 3

T = readtable(base_dir + "Hyperscanning_Log.csv", opts);  

T3 = T(1:length(dyad_lst), [1 2 4]);  % pilot number, dyad number, dyad type

% CHANGE. REMOVING 7TH ROW BC EXCL 7TH DYAD FOR NOW
T3([7],:) = [];

%% Load results
cd(deriv_dir)
load("IBS-values_excldyad7.mat")
load("IBS-values-block1_excldyad7.mat")
load("IBS-values-block2_excldyad7.mat")

% For dyads 3 and 4, the two PN runs only had 1 person guessing for each
% block. need to calc IBS using the first task blocks of the two runs 
PN_01_block1_avg_coh_dyad3 = IBS_block1.dyadic03.PN_01; % grab IBS results for run 1 PN task block 1 
PN_02_block1_avg_coh_dyad3 = IBS_block1.dyadic03.PN_02; % grab 2st task block for run 2
PN_01_block1_avg_coh_dyad4 = IBS_block1.dyadic04.PN_01; 
PN_02_block1_avg_coh_dyad4 = IBS_block1.dyadic04.PN_02;

for j = 1:10
    region_name = sprintf('Region%d',j);
    IBS_dyad3 = (PN_01_block1_avg_coh_dyad3.(region_name) + PN_02_block1_avg_coh_dyad3.(region_name)) / 2;
    IBS_dyad4 = (PN_01_block1_avg_coh_dyad4.(region_name) + PN_02_block1_avg_coh_dyad4.(region_name)) / 2;
    
    IBS_values.dyadic03.PN_01.(region_name) = IBS_dyad3;
    IBS_values.dyadic04.PN_01.(region_name) = IBS_dyad4;

end

%save('IBS-values_excl7_fixedcalc3and4.mat', 'IBS_values')

% IBS = (block1_avg_coh +block2_avg_coh)/2

%% Pull out PWA-SLP and PWA-NC into their own structs
results = IBS_values;
cfn_len = length(fieldnames(results)); % number of dyads
dyad_names = fieldnames(results);
IBS_PWA_SLP_tmp = {};  % patient clinician dyads
IBS_PWA_NC_tmp = {};  % patient healthy control dyads


for i =1:cfn_len
    curr_dyad = string(dyad_names(i))

    if contains(curr_dyad, string(T3.DyadNumber(i)))
        dyad_type = string(T3.DyadType(i));
        if contains(dyad_type, "PWA-SLP")
            IBS_PWA_SLP_tmp.(curr_dyad) = IBS_values.(curr_dyad);
        elseif contains(dyad_type, "PWA-NC")
            IBS_PWA_NC_tmp.(curr_dyad) = IBS_values.(curr_dyad);  % assign correct dyad data to new struct
        end
    end
end

%
% % manually input order
% Order them so that each matches the correct
desired_order_SLP = {'dyadic01','dyadic03','dyadic09','dyadic06','dyadic08'};
desired_order_NC  = {'dyadic02','dyadic04','dyadic05','dyadic10','dyadic11'}; 

IBS_PWA_SLP = orderfields(IBS_PWA_SLP_tmp, desired_order_SLP);
IBS_PWA_NC  = orderfields(IBS_PWA_NC_tmp,  desired_order_NC);


%% Look at overall average across all regions for each task 
% % Compare PWA-SLP to PWA-NC 
% task_lst = ["JS_01", "IR_01", "JR_01", "PN_01"];
% results = IBS_values;
% rows = [];  % to hold table rows
% 
% IBS_region_mean_SLP = {};
% IBS_region_sem_NC = {};
% for i = 1:length(dyad_lst)
%     curr_dyad = string(dyad_names(i));
%     for j = 1:length(task_lst)
%         curr_task = task_lst(j);
%         foo = IBS_values.(curr_dyad).(curr_task);
%         vals = struct2array(foo);
% 
%         foo_mean=mean(vals);
%         foo_sem = std(vals) / sqrt(numel(vals));
%         
%         IBS_region_mean.(curr_dyad).(curr_task) = foo_mean;
%         IBS_region_sem.(curr_dyad).(curr_task) = foo_sem;
%     end
% end
% 
%% Regions  (from MNI_corrected excel)

IFG = [1,6];
MFG=[2,7];
DLPFC = [3,8];
TPJ = [4,9];
SMG = [5,10];

LH =[1:5];
RH = [6:10];

%% Avg of all dyads for each ROI for each task

dyad_nms = string(fieldnames(IBS_values));

% Joint singing
[d_JS_pSLP, d_JS_pNC] = single_task_IBS_results('JS_01', IBS_PWA_SLP, IBS_PWA_NC, cfn_len, dyad_nms);  % reorganize data from struct for curr task
d_avg_hem_JS_pSLP = avg_hems(d_JS_pSLP); % average over the hemispheres
d_avg_hem_JS_pNC = avg_hems(d_JS_pNC); % rows are dyads, cols are ROI
[mean_JS_pSLP, sem_JS_pSLP] = group_mean_sem(d_avg_hem_JS_pSLP);  % calc group avg across dyads
[mean_JS_pNC, sem_JS_pNC] = group_mean_sem(d_avg_hem_JS_pNC);  % calc group avg across dyads
mean_JS = vertcat(mean_JS_pSLP, mean_JS_pNC); % concat into 1 matrix. top row = pSLP, bottom row = pNC
sem_JS = vertcat(sem_JS_pSLP, sem_JS_pNC);


% Ind reading
[d_IR_pSLP, d_IR_pNC] = single_task_IBS_results('IR_01', IBS_PWA_SLP, IBS_PWA_NC, cfn_len, dyad_nms);
% average over the hemispheres
d_avg_hem_IR_pSLP = avg_hems(d_IR_pSLP); % average over the hemispheres
d_avg_hem_IR_pNC = avg_hems(d_IR_pNC); % rows are dyads, cols are ROI
% calc group avg across dyads
[mean_IR_pSLP, sem_IR_pSLP] = group_mean_sem(d_avg_hem_IR_pSLP);  
[mean_IR_pNC, sem_IR_pNC] = group_mean_sem(d_avg_hem_IR_pNC);  % calc group avg across dyads
mean_IR = vertcat(mean_IR_pSLP, mean_IR_pNC); % concat into 1 matrix. top row = pSLP, bottom row = pNC
sem_IR = vertcat(sem_IR_pSLP, sem_IR_pNC);


%[mean_IR, sem_IR] = calc_group_mean_per_ROI_collapse_hems(d_IR_pSLP, d_IR_pNC); % calc mean across dyads for each ROI 

% Joint Reading
[d_JR_pSLP, d_JR_pNC] = single_task_IBS_results('JR_01', IBS_PWA_SLP, IBS_PWA_NC, cfn_len, dyad_nms);
d_avg_hem_JR_pSLP = avg_hems(d_JR_pSLP); % average over the hemispheres
d_avg_hem_JR_pNC = avg_hems(d_JR_pNC); % rows are dyads, cols are ROI
[mean_JR_pSLP, sem_JR_pSLP] = group_mean_sem(d_avg_hem_JR_pSLP);  % calc group avg across dyads
[mean_JR_pNC, sem_JR_pNC] = group_mean_sem(d_avg_hem_JR_pNC);  % calc group avg across dyads
mean_JR = vertcat(mean_JR_pSLP, mean_JR_pNC); % concat into 1 matrix. top row = pSLP, bottom row = pNC
sem_JR = vertcat(sem_JR_pSLP, sem_JR_pNC);
%[mean_JR, sem_JR] = calc_group_mean_per_ROI_collapse_hems(d_JR_pSLP, d_JR_pNC); % calc mean across dyads for each ROI 

% Pic Naming
[d_PN_pSLP, d_PN_pNC] = single_task_IBS_results('PN_01', IBS_PWA_SLP, IBS_PWA_NC, cfn_len, dyad_nms);
d_avg_hem_PN_pSLP = avg_hems(d_PN_pSLP); % average over the hemispheres
d_avg_hem_PN_pNC = avg_hems(d_PN_pNC); % rows are dyads, cols are ROI

[mean_PN_pSLP, sem_PN_pSLP] = group_mean_sem(d_avg_hem_PN_pSLP);  % calc group avg across dyads
[mean_PN_pNC, sem_PN_pNC] = group_mean_sem(d_avg_hem_PN_pNC);  % calc group avg across dyads

mean_PN = vertcat(mean_PN_pSLP, mean_PN_pNC); % concat into 1 matrix. top row = pSLP, bottom row = pNC
sem_PN = vertcat(sem_PN_pSLP, sem_PN_pNC);

disp('Group average across dyads calculated for each task and ROI')

%% Plot group avg over both hems on plot
%
% mean_task: top row = pSLP, bottom row = pNC
% cols in order: 
col_rois = {'IFG', 'MFG', 'DLPFC', 'TPJ', 'SMG'};
groups = {'patient-clinician', 'patient-control'}; % pSLP, pNC

task_lst_nms = ["Joint Singing", "Independent Reading", "Joint Reading", "Picture Naming"];

for i = 1:length(task_lst)
    figure()
    mean_task = mean_task_lst{i};
    sem_task = sem_task_lst{i};

    % Arrange as 5x2 so bar() makes 5 groups (ROIs) x 2 series (pSLP,nSLP)
    M = mean_task.';   % 5x2: columns = [pSLP nSLP]
    S = sem_task.';    % 5x2
    
    b = bar(M, 'grouped');    % grouped bars per ROI
    hold on
    
    % Add SEM error bars
    for k = 1:size(M,2)       % over the 2 series
        x = b(k).XEndPoints;  % x-locations of the bars for series k
        errorbar(x, M(:,k), S(:,k), 'k', 'linestyle', 'none', 'LineWidth', 1.2);
    end
    
    % Cosmetics
    ylim([0 0.45])
    set(gca, 'XTickLabel', col_rois, 'FontSize', 12);
    xlabel('ROI');
    ylabel('IBS (mean \pm SEM)');
    legend(groups, 'Location', 'northwest', 'Box', 'off');
    title(task_lst_nms(i) + "(mean \pm SEM)");
    grid off; box off;
    hold off
end

% Perform Wilcoxon sign ranked test

slp_keys = {'dyadic01','dyadic03','dyadic09','dyadic06','dyadic08'};
nc_keys  = {'dyadic02','dyadic04','dyadic05','dyadic10','dyadic11'}; 

task_fields = string(fieldnames(IBS_PWA_SLP.dyadic01));
roi_fields = string(fieldnames(IBS_PWA_SLP.dyadic01.PN_01));


stats = run_paired_tests(IBS_PWA_SLP, IBS_PWA_NC, slp_keys, nc_keys, task_fields, roi_fields);


%% take average of dyads over all rois

[avg_JS_pSLP2, sem_JS_pSLP2] = group_mean_sem_hems_sep(d_JS_pSLP);
[avg_IR_pSLP2, sem_IR_pSLP2]  = group_mean_sem_hems_sep(d_IR_pSLP);
[avg_JR_pSLP2, sem_JR_pSLP2] = group_mean_sem_hems_sep(d_JR_pSLP);
[avg_PN_pSLP2,sem_PN_pSLP2]  = group_mean_sem_hems_sep(d_PN_pSLP);

[avg_JS_pNC2, sem_JS_pNC2] = group_mean_sem_hems_sep(d_JS_pNC);
[avg_IR_pNC2, sem_IR_pNC2] = group_mean_sem_hems_sep(d_IR_pNC);
[avg_JR_pNC2, sem_JR_pNC2] = group_mean_sem_hems_sep(d_JR_pNC);
[avg_PN_pNC2, sem_PN_pNC2] = group_mean_sem_hems_sep(d_PN_pNC);

% Calc p value between hems
slp_keys = {'dyadic01','dyadic03','dyadic09','dyadic06','dyadic08'};
nc_keys  = {'dyadic02','dyadic04','dyadic05','dyadic10','dyadic11'};

task_fields = {'JS_01','IR_01','JR_01','PN_01'};

hemi_stats_slp = compare_hemispheres(IBS_PWA_SLP, slp_keys, task_fields); % taskXROI
hemi_stats_nc  = compare_hemispheres(IBS_PWA_NC,  nc_keys,  task_fields);

% grab p vals
p_JS_slp = hemi_stats_slp.pvals(1,:);  % 1x5
p_JS_nc  = hemi_stats_nc.pvals(1,:);   % 1x5
p_IR_slp = hemi_stats_slp.pvals(2,:);  % 1x5
p_IR_nc  = hemi_stats_nc.pvals(2,:);   % 1x5
p_JR_slp = hemi_stats_slp.pvals(3,:);  % 1x5
p_JR_nc  = hemi_stats_nc.pvals(3,:);   % 1x5
p_PN_slp = hemi_stats_slp.pvals(4,:);  % 1x5
p_PN_nc  = hemi_stats_nc.pvals(4,:);   % 1x5

% Plot left vs right

lr_save_path = save_dir + "leftvsright/";

cats=categorical(["Left-Brain"; "Right-Brain"]);

% Plot each task for patient-clinician
plot_leftvsright(avg_JS_pSLP2, sem_JS_pSLP2, cats, "PWA-clinician - Joint Singing", p_JS_slp, lr_save_path, "LvsR_slp_JS")
plot_leftvsright(avg_IR_pSLP2, sem_IR_pSLP2, cats, 'PWA-clinician - Independent Reading', p_IR_slp, lr_save_path, "LvsR_slp_IR")
plot_leftvsright(avg_JR_pSLP2, sem_JR_pSLP2, cats, 'PWA-clinician - Joint Reading',p_JR_slp, lr_save_path, "LvsR_slp_JR")
plot_leftvsright(avg_PN_pSLP2, sem_PN_pSLP2, cats, 'PWA-clinician - Picture Guessing',p_PN_slp, lr_save_path, "LvsR_slp_PN")

% Plot each task for patient-stranger
plot_leftvsright(avg_JS_pNC2, sem_JS_pNC2, cats, "PWA-stranger - Joint Singing",p_JS_nc, lr_save_path, "LvsR_nc_JS")
plot_leftvsright(avg_IR_pNC2, sem_IR_pNC2, cats, 'PWA-stranger - Independent Reading',p_IR_nc, lr_save_path, "LvsR_nc_IR")
plot_leftvsright(avg_JR_pNC2, sem_JR_pNC2, cats, 'PWA-stranger - Joint Reading',p_JR_nc, lr_save_path, "LvsR_nc_JR")
plot_leftvsright(avg_PN_pNC2, sem_PN_pNC2, cats, 'PWA-stranger - Picture Guessing',p_PN_nc, lr_save_path, "LvsR_nc_PN")


%% Global comparisons -- collapse across all left ROIs and all right ROIs
d_avg_roiperhem_JS_pSLP = avg_rois_per_hem(d_JS_pSLP); % dyad X hemisphere
d_avg_roiperhem_IR_pSLP = avg_rois_per_hem(d_IR_pSLP);
d_avg_roiperhem_JR_pSLP = avg_rois_per_hem(d_JR_pSLP);
d_avg_roiperhem_PN_pSLP = avg_rois_per_hem(d_PN_pSLP);

d_avg_roiperhem_JS_pNC = avg_rois_per_hem(d_JS_pNC);
d_avg_roiperhem_IR_pNC = avg_rois_per_hem(d_IR_pNC);
d_avg_roiperhem_JR_pNC = avg_rois_per_hem(d_JR_pNC);
d_avg_roiperhem_PN_pNC = avg_rois_per_hem(d_PN_pNC);

%
% calc group average, sem, p vals and plot
%
lr_save_path = save_dir + "global_hems/";
plot_hemi_group(d_avg_roiperhem_JS_pSLP, "PWA-clinician - Joint Singing", lr_save_path, "LvsR_slp_JS");
plot_hemi_group(d_avg_roiperhem_IR_pSLP, "PWA-clinician - Independent Reading", lr_save_path, "LvsR_slp_IR");
plot_hemi_group(d_avg_roiperhem_JR_pSLP, "PWA-clinician - Joint Reading", lr_save_path, "LvsR_slp_JR");
plot_hemi_group(d_avg_roiperhem_PN_pSLP, "PWA-clinician - Picture Guessing",lr_save_path, "LvsR_slp_PN");

plot_hemi_group(d_avg_roiperhem_JS_pNC, "PWA-stranger - Joint Singing", lr_save_path, "LvsR_nc_JS");
plot_hemi_group(d_avg_roiperhem_IR_pNC, "PWA-stranger - Independent Reading", lr_save_path, "LvsR_nc_IR");
plot_hemi_group(d_avg_roiperhem_JR_pNC, "PWA-stranger - Joint Reading", lr_save_path, "LvsR_nc_JR");
plot_hemi_group(d_avg_roiperhem_PN_pNC, "PWA-stranger - Picture Guessing",lr_save_path, "LvsR_nc_PN");


% 
% % get p vals
% stats_JS_slp = paired_hemi_stats(d_avg_roiperhem_JS_pSLP); p_JS_slp = stats_JS_slp.p;
% stats_IR_slp = paired_hemi_stats(d_avg_roiperhem_IR_pSLP); p_IR_slp = stats_IR_slp.p;
% stats_JR_slp = paired_hemi_stats(d_avg_roiperhem_JR_pSLP); p_JR_slp = stats_JR_slp.p;
% stats_PN_slp = paired_hemi_stats(d_avg_roiperhem_PN_pSLP); p_PN_slp = stats_PN_slp.p;
% 
% stats_JS_nc = paired_hemi_stats(d_avg_roiperhem_JS_pNC); p_JS_nc = stats_JS_nc.p;
% stats_IR_nc = paired_hemi_stats(d_avg_roiperhem_IR_pNC); p_IR_nc = stats_IR_nc.p;
% stats_JR_nc = paired_hemi_stats(d_avg_roiperhem_JR_pNC); p_JR_nc = stats_JR_nc.p;
% stats_PN_nc = paired_hemi_stats(d_avg_roiperhem_PN_pNC); p_PN_nc = stats_PN_nc.p;
% 
% % get group mean
% [mean_LR_JS_pSLP, sem_LR_JS_pSLP] = group_mean_sem(d_avg_roiperhem_JS_pSLP);
% [mean_LR_IR_pSLP, sem_LR_IR_pSLP] = group_mean_sem(d_avg_roiperhem_IR_pSLP);
% [mean_LR_JR_pSLP, sem_LR_JR_pSLP] = group_mean_sem(d_avg_roiperhem_JR_pSLP);
% [mean_LR_PN_pSLP, sem_LR_PN_pSLP] = group_mean_sem(d_avg_roiperhem_PN_pSLP);
% 
% [mean_LR_JS_pNC, sem_LR_JS_pNC] = group_mean_sem(d_avg_roiperhem_JS_pNC);
% [mean_LR_IR_pNC, sem_LR_IR_pNC] = group_mean_sem(d_avg_roiperhem_IR_pNC);
% [mean_LR_JR_pNC, sem_LR_JR_pNC] = group_mean_sem(d_avg_roiperhem_JR_pNC);
% [mean_LR_PN_pNC, sem_LR_PN_pNC] = group_mean_sem(d_avg_roiperhem_PN_pNC);
% 

%% Now compare across tasks for each hemisphere

tasks = {'JS','IR','JR','PN'};
tasks = {'Joint Singing','Independent Reading','Joint Reading','Picture Guessing'};


% Pack SLP & NC mats for each task
mats_slp = { d_avg_roiperhem_JS_pSLP, d_avg_roiperhem_IR_pSLP, ...  cells of dyad X hemisphere for each task
             d_avg_roiperhem_JR_pSLP, d_avg_roiperhem_PN_pSLP };

mats_nc  = { d_avg_roiperhem_JS_pNC,  d_avg_roiperhem_IR_pNC,  ...
             d_avg_roiperhem_JR_pNC,  d_avg_roiperhem_PN_pNC  };

% LEFT hemisphere across tasks: SLP vs NC per task
plot_SLPvsNC_across_tasks_left_or_right(mats_slp, mats_nc, tasks, 1, 'Left Hemisphere across tasks');
% RIGHT hemisphere across tasks: SLP vs NC per task
plot_SLPvsNC_across_tasks_left_or_right(mats_slp, mats_nc, tasks, 2, 'Right Hemisphere across tasks');

%% Compare across tasks for EACH DYAD

slp_keys = {'dyadic01','dyadic03','dyadic09','dyadic06','dyadic08'};
nc_keys  = {'dyadic02','dyadic04','dyadic05','dyadic10','dyadic11'};
tasks = {'Joint Singing','Independent Reading','Joint Reading','Picture Guessing'};


[d_avg_roiperhem_JS_pSLP, d_sem_roiperhem_JS_pSLP] = avg_sem_rois_per_hem(d_JS_pSLP); % dyad X hemisphere
[d_avg_roiperhem_IR_pSLP, d_sem_roiperhem_IR_pSLP] = avg_sem_rois_per_hem(d_IR_pSLP);
[d_avg_roiperhem_JR_pSLP, d_sem_roiperhem_JR_pSLP] = avg_sem_rois_per_hem(d_JR_pSLP);
[d_avg_roiperhem_PN_pSLP, d_sem_roiperhem_PN_pSLP] = avg_sem_rois_per_hem(d_PN_pSLP);

[d_avg_roiperhem_JS_pNC, d_sem_roiperhem_JS_pNC] = avg_sem_rois_per_hem(d_JS_pNC);
[d_avg_roiperhem_IR_pNC, d_sem_roiperhem_IR_pNC] = avg_sem_rois_per_hem(d_IR_pNC);
[d_avg_roiperhem_JR_pNC, d_sem_roiperhem_JR_pNC] = avg_sem_rois_per_hem(d_JR_pNC);
[d_avg_roiperhem_PN_pNC, d_sem_roiperhem_PN_pNC] = avg_sem_rois_per_hem(d_PN_pNC);

% plot each dyad pair (SLP-NC) on same plot with all tasks on x axis
% do it for each hemisphere (separate plots for dyad pairs and hemispheres)
slp_avg = {d_avg_roiperhem_JS_pSLP, d_avg_roiperhem_IR_pSLP, ...
           d_avg_roiperhem_JR_pSLP, d_avg_roiperhem_PN_pSLP};
nc_avg  = {d_avg_roiperhem_JS_pNC,  d_avg_roiperhem_IR_pNC, ...
           d_avg_roiperhem_JR_pNC,  d_avg_roiperhem_PN_pNC};

slp_sem = {d_sem_roiperhem_JS_pSLP, d_sem_roiperhem_IR_pSLP, ...
           d_sem_roiperhem_JR_pSLP, d_sem_roiperhem_PN_pSLP};
nc_sem  = {d_sem_roiperhem_JS_pNC,  d_sem_roiperhem_IR_pNC, ...
           d_sem_roiperhem_JR_pNC,  d_sem_roiperhem_PN_pNC};

% 
% compute ttest
%
hems  = {'Left','Right'};
nTasks = numel(tasks);
nHems  = size(slp_avg{1},2);

% Store results
pvals_slp = cell(nHems,1);
pvals_nc  = cell(nHems,1);

for h = 1:nHems
    % Build dyad × task matrices
    X_slp = nan(size(slp_avg{1},1), nTasks);
    X_nc  = nan(size(nc_avg{1},1), nTasks);
    for t = 1:nTasks
        X_slp(:,t) = slp_avg{t}(:,h);  % dyads × tasks
        X_nc(:,t)  = nc_avg{t}(:,h);
    end

    % Paired t-tests across dyads (task vs task, within group)
    fprintf('\n=== %s Hemisphere ===\n', hems{h});
    fprintf('--- SLP group ---\n');
    for i = 1:nTasks
        for j = i+1:nTasks
            [~,p,~,stats] = ttest(X_slp(:,i), X_slp(:,j));
            fprintf('%s vs %s: t(%d)=%.2f, p=%.4f\n', ...
                tasks{i}, tasks{j}, stats.df, stats.tstat, p);
            pvals_slp{h}(i,j) = p;
        end
    end

    fprintf('--- NC group ---\n');
    for i = 1:nTasks
        for j = i+1:nTasks
            [~,p,~,stats] = ttest(X_nc(:,i), X_nc(:,j));
            fprintf('%s vs %s: t(%d)=%.2f, p=%.4f\n', ...
                tasks{i}, tasks{j}, stats.df, stats.tstat, p);
            pvals_nc{h}(i,j) = p;
        end
    end
end


nTasks = numel(tasks);
nHems  = numel(pvals_nc);   % usually 2

for h = 1:nHems
    % --- NC ---
    U = pvals_nc{h};                      % expected size: (nTasks-1) x nTasks (upper-triangular block)
    P = nan(nTasks);                      % final nTasks x nTasks
    % fill upper triangle from rectangular U
    for i = 1:min(nTasks-1, size(U,1))
        for j = i+1:min(nTasks, size(U,2))
            P(i,j) = U(i,j);
        end
    end
    % mirror to lower triangle; keep diagonal NaN
    P = P + P.';
    P(1:nTasks+1:end) = NaN;
    pvals_nc{h} = P;

    % --- SLP ---
    U = pvals_slp{h};
    P = nan(nTasks);
    for i = 1:min(nTasks-1, size(U,1))
        for j = i+1:min(nTasks, size(U,2))
            P(i,j) = U(i,j);
        end
    end
    P = P + P.';
    P(1:nTasks+1:end) = NaN;
    pvals_slp{h} = P;
end


%

%tasks = {'JS','IR','JR','PN'};
hems  = {'Left','Right'};
pilot_labels = ['1', '2', '3', '4', '6'];

plot_bar_eachDyad(slp_avg, nc_avg, slp_sem, nc_sem, tasks, hems, pilot_labels, pvals_slp, pvals_nc);


%%
% 
% plot_SLPvsNC_across_tasks_left_or_right(mats_slp, mats_nc, task_labels, hemi_col, title_nm)
% 
% % Left hemisphere dyad-level per task (SLP)
% L_JS_slp = d_avg_roiperhem_JS_pSLP(:,1);
% L_IR_slp = d_avg_roiperhem_IR_pSLP(:,1);
% L_JR_slp = d_avg_roiperhem_JR_pSLP(:,1);
% L_PN_slp = d_avg_roiperhem_PN_pSLP(:,1);
% 
% L_JS_nc = d_avg_roiperhem_JS_pNC(:,1);
% L_IR_nc = d_avg_roiperhem_IR_pNC(:,1);
% L_JR_nc = d_avg_roiperhem_JR_pNC(:,1);
% L_PN_nc = d_avg_roiperhem_PN_pNC(:,1);
% 
% % Right hemisphere dyad-level per task (SLP)
% R_JS_slp = d_avg_roiperhem_JS_pSLP(:,2);
% R_IR_slp = d_avg_roiperhem_IR_pSLP(:,2);
% R_JR_slp = d_avg_roiperhem_JR_pSLP(:,2);
% R_PN_slp = d_avg_roiperhem_PN_pSLP(:,2);
% 
% R_JS_nc = d_avg_roiperhem_JS_pNC(:,2);
% R_IR_nc = d_avg_roiperhem_IR_pNC(:,2);
% R_JR_nc = d_avg_roiperhem_JR_pNC(:,2);
% R_PN_nc = d_avg_roiperhem_PN_pNC(:,2);


%% Compare tasks across ROIs (but separate for left and right ROIs)
% ROIs in order:
% IFG, MFG, DLPFC, TPJ, SMG

ROI_nms = {'IFG', 'MFG', 'DLPFC', 'TPJ', 'SMG'};
save_path = fullfile(save_dir, "compare_tasks");

% Your existing inputs
slp_keys = {'dyadic01','dyadic03','dyadic09','dyadic06','dyadic08'};
nc_keys  = {'dyadic02','dyadic04','dyadic05','dyadic10','dyadic11'};
task_fields = {'JS_01','IR_01','JR_01','PN_01'};            % use your exact task field names
%task_labels = {'JS','IR','JR','PN'};
task_labels = {'Joint Singing','Independent Reading','Joint Reading','Picture Guessing'};

% %
% % IFG Left hemisphere
% roi_idx = 1; hemi = 'L';
% % Build dyad×task matrices
% X_slp = roi_task_matrix(IBS_PWA_SLP, slp_keys, task_fields, roi_idx, hemi);  % [nDyads × 4 tasks]
% X_nc  = roi_task_matrix(IBS_PWA_NC,  nc_keys,  task_fields, roi_idx, hemi);  % [nDyads × 4] tasks]
% % Compute RM stats per dyad type
% stats_slp = rm_task_tests(X_slp, task_labels);
% stats_nc  = rm_task_tests(X_nc,  task_labels);
% % Plot bars (means±SEM) + print significant pairs (FDR) for each group
% plot_tasks_for_roi(stats_slp.mu, stats_slp.sem, task_labels, 'PWA-clinician - Left IFG', stats_slp.p_pair, roi_idx, save_path, ROI_nms);
% plot_tasks_for_roi(stats_nc.mu,  stats_nc.sem,  task_labels, 'PWA-stranger - Left IFG', stats_nc.p_pair, roi_idx, save_path, ROI_nms);
% 
% %
% % MFG Left hemisphere
% roi_idx = 2; hemi = 'L';
% % Build dyad×task matrices
% X_slp = roi_task_matrix(IBS_PWA_SLP, slp_keys, task_fields, roi_idx, hemi);  % [nDyads × 4 tasks]
% X_nc  = roi_task_matrix(IBS_PWA_NC,  nc_keys,  task_fields, roi_idx, hemi);  % [nDyads × 4] tasks]
% % Compute RM stats per dyad type
% stats_slp = rm_task_tests(X_slp, task_labels);
% stats_nc  = rm_task_tests(X_nc,  task_labels);
% % Plot bars (means±SEM) + print significant pairs (FDR) for each group
% plot_tasks_for_roi(stats_slp.mu, stats_slp.sem, task_labels, 'PWA-clinician - Left MFG', stats_slp.p_pair, roi_idx, save_path, ROI_nms);
% plot_tasks_for_roi(stats_nc.mu,  stats_nc.sem,  task_labels, 'PWA-stranger - Left MFG', stats_nc.p_pair, roi_idx, save_path, ROI_nms);
% 
% %
% % DLPFC Left hemisphere
% roi_idx = 3; hemi = 'L';
% % Build dyad×task matrices
% X_slp = roi_task_matrix(IBS_PWA_SLP, slp_keys, task_fields, roi_idx, hemi);  % [nDyads × 4 tasks]
% X_nc  = roi_task_matrix(IBS_PWA_NC,  nc_keys,  task_fields, roi_idx, hemi);  % [nDyads × 4] tasks]
% % Compute RM stats per dyad type
% stats_slp = rm_task_tests(X_slp, task_labels);
% stats_nc  = rm_task_tests(X_nc,  task_labels);
% % Plot bars (means±SEM) + print significant pairs (FDR) for each group
% plot_tasks_for_roi(stats_slp.mu, stats_slp.sem, task_labels, 'PWA-clinician - Left DLPFC', stats_slp.p_pair, roi_idx, save_path, ROI_nms);
% plot_tasks_for_roi(stats_nc.mu,  stats_nc.sem,  task_labels, 'PWA-stranger - Left DLPFC', stats_nc.p_pair, roi_idx, save_path, ROI_nms);
% 
% %
% % TPJ Left hemisphere
% roi_idx = 4; hemi = 'L';
% % Build dyad×task matrices
% X_slp = roi_task_matrix(IBS_PWA_SLP, slp_keys, task_fields, roi_idx, hemi);  % [nDyads × 4 tasks]
% X_nc  = roi_task_matrix(IBS_PWA_NC,  nc_keys,  task_fields, roi_idx, hemi);  % [nDyads × 4] tasks]
% % Compute RM stats per dyad type
% stats_slp = rm_task_tests(X_slp, task_labels);
% stats_nc  = rm_task_tests(X_nc,  task_labels);
% % Plot bars (means±SEM) + print significant pairs (FDR) for each group
% plot_tasks_for_roi(stats_slp.mu, stats_slp.sem, task_labels, 'PWA-clinician - Left TPJ', stats_slp.p_pair, roi_idx, save_path, ROI_nms);
% plot_tasks_for_roi(stats_nc.mu,  stats_nc.sem,  task_labels, 'PWA-stranger - Left TPJ', stats_nc.p_pair, roi_idx, save_path, ROI_nms);
% 
% %
% % SMG Left hemisphere
% roi_idx = 5; hemi = 'L';
% % Build dyad×task matrices
% X_slp = roi_task_matrix(IBS_PWA_SLP, slp_keys, task_fields, roi_idx, hemi);  % [nDyads × 4 tasks]
% X_nc  = roi_task_matrix(IBS_PWA_NC,  nc_keys,  task_fields, roi_idx, hemi);  % [nDyads × 4] tasks]
% % Compute RM stats per dyad type
% stats_slp = rm_task_tests(X_slp, task_labels);
% stats_nc  = rm_task_tests(X_nc,  task_labels);
% % Plot bars (means±SEM) + print significant pairs (FDR) for each group
% plot_tasks_for_roi(stats_slp.mu, stats_slp.sem, task_labels, 'PWA-clinician - Left SMG', stats_slp.p_pair, roi_idx, save_path, ROI_nms);
% plot_tasks_for_roi(stats_nc.mu,  stats_nc.sem,  task_labels, 'PWA-stranger - Left SMG', stats_nc.p_pair, roi_idx, save_path, ROI_nms);
% 
% %
% % RIGHT hem 
% %
% %
% % IFG RH
% roi_idx = 1; hemi = 'R';
% % Build dyad×task matrices
% X_slp = roi_task_matrix(IBS_PWA_SLP, slp_keys, task_fields, roi_idx, hemi);  % [nDyads × 4 tasks]
% X_nc  = roi_task_matrix(IBS_PWA_NC,  nc_keys,  task_fields, roi_idx, hemi);  % [nDyads × 4] tasks]
% % Compute RM stats per dyad type
% stats_slp = rm_task_tests(X_slp, task_labels);
% stats_nc  = rm_task_tests(X_nc,  task_labels);
% % Plot bars (means±SEM) + print significant pairs (FDR) for each group
% plot_tasks_for_roi(stats_slp.mu, stats_slp.sem, task_labels, 'PWA-clinician - Right IFG', stats_slp.p_pair, roi_idx, save_path, ROI_nms);
% plot_tasks_for_roi(stats_nc.mu,  stats_nc.sem,  task_labels, 'PWA-stranger - Right IFG', stats_nc.p_pair, roi_idx, save_path, ROI_nms);
% 
% %
% % MFG Left hemisphere
% roi_idx = 2; hemi = 'R';
% % Build dyad×task matrices
% X_slp = roi_task_matrix(IBS_PWA_SLP, slp_keys, task_fields, roi_idx, hemi);  % [nDyads × 4 tasks]
% X_nc  = roi_task_matrix(IBS_PWA_NC,  nc_keys,  task_fields, roi_idx, hemi);  % [nDyads × 4] tasks]
% % Compute RM stats per dyad type
% stats_slp = rm_task_tests(X_slp, task_labels);
% stats_nc  = rm_task_tests(X_nc,  task_labels);
% % Plot bars (means±SEM) + print significant pairs (FDR) for each group
% plot_tasks_for_roi(stats_slp.mu, stats_slp.sem, task_labels, 'PWA-clinician - Right MFG', stats_slp.p_pair, roi_idx, save_path, ROI_nms);
% plot_tasks_for_roi(stats_nc.mu,  stats_nc.sem,  task_labels, 'PWA-stranger - Right MFG', stats_nc.p_pair, roi_idx, save_path, ROI_nms);
% 
% %
% % DLPFC Left hemisphere
% roi_idx = 3; hemi = 'R';
% % Build dyad×task matrices
% X_slp = roi_task_matrix(IBS_PWA_SLP, slp_keys, task_fields, roi_idx, hemi);  % [nDyads × 4 tasks]
% X_nc  = roi_task_matrix(IBS_PWA_NC,  nc_keys,  task_fields, roi_idx, hemi);  % [nDyads × 4] tasks]
% % Compute RM stats per dyad type
% stats_slp = rm_task_tests(X_slp, task_labels);
% stats_nc  = rm_task_tests(X_nc,  task_labels);
% % Plot bars (means±SEM) + print significant pairs (FDR) for each group
% plot_tasks_for_roi(stats_slp.mu, stats_slp.sem, task_labels, 'PWA-clinician - Right DLPFC', stats_slp.p_pair, roi_idx, save_path, ROI_nms);
% plot_tasks_for_roi(stats_nc.mu,  stats_nc.sem,  task_labels, 'PWA-stranger - Right DLPFC', stats_nc.p_pair, roi_idx, save_path, ROI_nms);
% 
% %
% % TPJ Left hemisphere
% roi_idx = 4; hemi = 'R';
% % Build dyad×task matrices
% X_slp = roi_task_matrix(IBS_PWA_SLP, slp_keys, task_fields, roi_idx, hemi);  % [nDyads × 4 tasks]
% X_nc  = roi_task_matrix(IBS_PWA_NC,  nc_keys,  task_fields, roi_idx, hemi);  % [nDyads × 4] tasks]
% % Compute RM stats per dyad type
% stats_slp = rm_task_tests(X_slp, task_labels);
% stats_nc  = rm_task_tests(X_nc,  task_labels);
% % Plot bars (means±SEM) + print significant pairs (FDR) for each group
% plot_tasks_for_roi(stats_slp.mu, stats_slp.sem, task_labels, 'PWA-clinician - Right TPJ', stats_slp.p_pair, roi_idx, save_path, ROI_nms);
% plot_tasks_for_roi(stats_nc.mu,  stats_nc.sem,  task_labels, 'PWA-stranger - Right TPJ', stats_nc.p_pair, roi_idx, save_path, ROI_nms);
% 
% %
% % SMG Left hemisphere
% roi_idx = 5; hemi = 'R';
% % Build dyad×task matrices
% X_slp = roi_task_matrix(IBS_PWA_SLP, slp_keys, task_fields, roi_idx, hemi);  % [nDyads × 4 tasks]
% X_nc  = roi_task_matrix(IBS_PWA_NC,  nc_keys,  task_fields, roi_idx, hemi);  % [nDyads × 4] tasks]
% % Compute RM stats per dyad type
% stats_slp = rm_task_tests(X_slp, task_labels);
% stats_nc  = rm_task_tests(X_nc,  task_labels);
% % Plot bars (means±SEM) + print significant pairs (FDR) for each group
% plot_tasks_for_roi(stats_slp.mu, stats_slp.sem, task_labels, 'PWA-clinician - Right SMG', stats_slp.p_pair, roi_idx, save_path, ROI_nms);
% plot_tasks_for_roi(stats_nc.mu,  stats_nc.sem,  task_labels, 'PWA-stranger - Right SMG', stats_nc.p_pair, roi_idx, save_path, ROI_nms);
% 
% FolderName = fullfile(save_path, "png");   % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   %FigName   = get(FigHandle, 'Name');
%   FigName   = ['Fig' num2str(iFig)];
%   %savefig(FigHandle, fullfile(FolderName, FigName, '.png'));
%   saveas(FigHandle,fullfile(FolderName, [FigName '.png'])); %Specify format for the figure
% end

%% diff fig format - Compare tasks across ROIs (but separate for left and right ROIs)
% ROI labels (optional for figure titles)
roi_labels = {'IFG','MFG','DLPFC','TPJ','SMG'};
hemi_labels = struct('L','Left','R','Right');

save_path = fullfile(save_dir, "compare_tasks_diff_format");

% Task names
%task_fields = {'JS','IR','JR','PN'};
%task_labels = {'JS','IR','JR','PN'};

for roi_idx = 1:5                      % loop ROIs
    for hemi_char = ['L','R']          % loop hemispheres
        hemi = char(hemi_char);        % 'L' or 'R'
        
        % Build dyad×task matrices for SLP and NC
        X_slp = roi_task_matrix(IBS_PWA_SLP, slp_keys, task_fields, roi_idx, hemi);
        X_nc  = roi_task_matrix(IBS_PWA_NC,  nc_keys,  task_fields, roi_idx, hemi);

        % Compute task-wise stats (means, SEMs, pairwise tests)
        stats_slp = rm_task_tests(X_slp, task_labels);
        stats_nc  = rm_task_tests(X_nc,  task_labels);

        % Figure title
        fig_title = sprintf('%s %s', hemi_labels.(hemi), roi_labels{roi_idx});

        % Plot both groups with task-comparison stars for BOTH groups
        plot_tasks_SLP_NC_with_two_group_sigs(stats_slp, stats_nc, task_labels, ...
            fig_title, 'UseFDR', false, 'OnlyAdjacent', false, 'MaxPairsEach', 6);

        % (Optional) save figures
        save_root = fullfile(save_dir, 'taskwise');
        mkdir(save_root);
        fname = sprintf('%s_%s_taskwise.png', hemi, roi_labels{roi_idx});
        exportgraphics(gcf, fullfile(save_root, fname), 'Resolution', 300);
    end
end

%% Compare SLP-NC dyads w all rois incl left and right
slp_keys = {'dyadic01','dyadic03','dyadic09','dyadic06','dyadic08'};
nc_keys  = {'dyadic02','dyadic04','dyadic05','dyadic10','dyadic11'};
task_fields = {'JS_01','IR_01','JR_01','PN_01'};
rois = {'IFG','MFG','DLPFC','TPJ','SMG'};

X_slp_JS = extract_task_matrix(IBS_PWA_SLP, slp_keys, 'JS_01');  % [nDyads×10]
X_nc_JS  = extract_task_matrix(IBS_PWA_NC,  nc_keys,  'JS_01');
task_label = 'Joint Singing';
for r = 1:numel(rois)
    plot_SLPvsNC_byROI_hemi(X_slp_JS, X_nc_JS, rois{r}, task_label);
end

X_slp_IR = extract_task_matrix(IBS_PWA_SLP, slp_keys, 'IR_01');  % [nDyads×10]
X_nc_IR  = extract_task_matrix(IBS_PWA_NC,  nc_keys,  'IR_01');
task_label = 'Independent Reading';
for r = 1:numel(rois)
    plot_SLPvsNC_byROI_hemi(X_slp_IR, X_nc_IR, rois{r}, task_label);
end


X_slp_JR = extract_task_matrix(IBS_PWA_SLP, slp_keys, 'JR_01');  % [nDyads×10]
X_nc_JR  = extract_task_matrix(IBS_PWA_NC,  nc_keys,  'JR_01');
task_label = 'Joint Reading';
for r = 1:numel(rois)
    plot_SLPvsNC_byROI_hemi(X_slp_JR, X_nc_JR, rois{r}, task_label);
end


X_slp_PN = extract_task_matrix(IBS_PWA_SLP, slp_keys, 'PN_01');  % [nDyads×10]
X_nc_PN  = extract_task_matrix(IBS_PWA_NC,  nc_keys,  'PN_01');
task_label = 'Picture Guessing';
for r = 1:numel(rois)
    plot_SLPvsNC_byROI_hemi(X_slp_PN, X_nc_PN, rois{r}, task_label);
end

%%
save_path = fullfile(save_dir, "individ_dyads_LR");

FolderName = save_path;   % Your destination folder
if ~isfolder(FolderName) % Check if the folder does not exist
    mkdir(FolderName);   % Create the folder
end
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  %FigName   = get(FigHandle, 'Name');
  FigName   = ['Fig' num2str(iFig)];
  %savefig(FigHandle, fullfile(FolderName, FigName, '.png'));
  saveas(FigHandle,fullfile(FolderName, [FigName '.fig'])); %Specify format for the figure
end



%%
% % take the avg of all ROI for a single task 
% for a = 1:5
%     for m = 1:cfn_len
%     meanvec(m) = d_JS{m}(1,a);
%     rJS_avg(1,a) = mean(meanvec, 'omitnan');
%     end
% end
% for a = 1:5
%     for m = 1:cfn_len
%     meanvec(m) = d_JS{m}(2,a);
%     rJS_avg(2,a) = mean(meanvec, 'omitnan');
%     end
% end

%% Stats
% 
% mean_task_lst = {mean_JS, mean_IR, mean_JR, mean_PN};
% sem_task_lst = {sem_JS, sem_IR, sem_JR, sem_PN};
% % Perform paired ttest --- not good bc we only have group avgs
% % put the 2 dyad types into 2 arrays
% nTasks = numel(mean_task_lst);
% nROI   = size(mean_task_lst{1}, 2);
% 
% % Preallocate
% SLP_vals     = nan(nTasks, nROI);   % top row
% Stranger_vals = nan(nTasks, nROI);  % bottom row
% 
% for t = 1:nTasks
%     tmp = mean_task_lst{t};
%     SLP_vals(t,:)      = tmp(1,:);  % row 1
%     Stranger_vals(t,:) = tmp(2,:);  % row 2
% end
% 
% % 
% % Paired ttests ROI by ROI
% nTasks = numel(mean_task_lst);
% nROI   = size(mean_task_lst{1}, 2);
% 
% pvals_hemavg_roi = nan(1, nROI);
% stats_hemavg_roi = cell(1, nROI);
% for r = 1:nROI
%     [~, p, ~, stat_tmp] = ttest(SLP_vals(:,r), Stranger_vals(:,r));
%     pvals_hemavg_roi(r) = p;
%     stats_hemavg_roi{r} = stat_tmp;  % contains t-statistic, df, etc.
% end
% 
% pvals_hemavg_task = nan(1, nTasks);
% stats_hemavg_task = cell(1, nTasks);
% for t = 1:nTasks
%     [~, p, ~, stat_tmp] = ttest(SLP_vals(t,:), Stranger_vals(t,:));
%     pvals_hemavg_task(t) = p;
%     stats_hemavg_task{t} = stat_tmp;  % contains t-statistic, df, etc.
% end

%%
% %% Perform paired ttest --- not good bc we only have group avgs
% % put the 2 dyad types into 2 arrays
% nTasks = numel(mean_task_lst);
% nROI   = size(mean_task_lst{1}, 2);
% 
% % Preallocate
% SLP_vals     = nan(nTasks, nROI);   % top row
% Stranger_vals = nan(nTasks, nROI);  % bottom row
% 
% for t = 1:nTasks
%     tmp = mean_task_lst{t};
%     SLP_vals(t,:)      = tmp(1,:);  % row 1
%     Stranger_vals(t,:) = tmp(2,:);  % row 2
% end
% 
% 
% % 
% % Paired ttests ROI by ROI
% nTasks = numel(mean_task_lst);
% nROI   = size(mean_task_lst{1}, 2);
% 
% % Preallocate
% SLP_vals     = nan(nTasks, nROI);   % top row
% Stranger_vals = nan(nTasks, nROI);  % bottom row
% 
% for t = 1:nTasks
%     tmp = mean_task_lst{t};
%     SLP_vals(t,:)      = tmp(1,:);  % row 1
%     Stranger_vals(t,:) = tmp(2,:);  % row 2
% end

%% plot each dyad pair

function plot_bar_eachDyad(d_avg_slp, d_avg_nc, d_sem_slp, d_sem_nc, task_labels, hem_labels, pilot_labels, pvals_slp, pvals_nc)
% d_avg_slp/nc: cell arrays {1×nTasks}, each element [nDyads × 2]
% d_sem_slp/nc: same structure, SEM across ROIs per dyad × hemisphere
% task_labels:  e.g. {'JS','IR','JR','PN'}
% hem_labels:   e.g. {'Left','Right'}

    nTasks = numel(task_labels);
    nDyads = size(d_avg_slp{1},1);  
    nHems  = size(d_avg_slp{1},2);

    for h = 1:nHems
        for d = 1:nDyads  % num pilots
            mu_slp  = nan(1,nTasks);
            mu_nc   = nan(1,nTasks);
            sem_slp = nan(1,nTasks);
            sem_nc  = nan(1,nTasks);

            for t = 1:nTasks
                mu_slp(t)  = d_avg_slp{t}(d,h);
                mu_nc(t)   = d_avg_nc{t}(d,h);
                sem_slp(t) = d_sem_slp{t}(d,h);
                sem_nc(t)  = d_sem_nc{t}(d,h);
            end
            
            % New figure per dyad × hemisphere
            figure('Name',['Pilot ' pilot_labels(d) ' - ' hem_labels{h}]);
            bh = bar([mu_slp; mu_nc]','grouped'); 
            hold on;

            % Error bars
            ngroups = nTasks; nbars = 2;
            groupwidth = min(0.8, nbars/(nbars+1.5));
            for i = 1:nbars
                x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
                if i == 1
                    errorbar(x, mu_slp, sem_slp, 'k', 'linestyle','none','LineWidth',1.2, 'HandleVisibility', 'off');
                else
                    errorbar(x, mu_nc, sem_nc, 'k', 'linestyle','none','LineWidth',1.2, 'HandleVisibility', 'off');
                end
            end

            set(gca,'XTick',1:nTasks,'XTickLabel',task_labels, 'FontSize', 12);
            title(['Pilot ' pilot_labels(d) ' - ' hem_labels{h} ' Hemisphere'], 'FontSize', 14);
            %ylabel('Mean ± SEM across ROIs');
            ylabel('Interbrain Synchrony', 'FontSize', 14);
            ylim([0 0.5]);
            legend({'PWA-clinician','PWA-stranger'},'Location', 'northeast', 'FontSize', 10);
        end

        % --- Plot stars from pvals ---
        % SLP = left bar (x = task - 0.15), NC = right bar (x = task + 0.15)
        offset = 0.15;
        for i = 1:nTasks
            for j = i+1:nTasks
                % SLP stars
                p = pvals_slp{h}(i,j);
                if ~isnan(p) && p < 0.05
                    y = max([mu_slp(i)+sem_slp_vec(i), mu_slp(j)+sem_slp_vec(j)]) * 1.2;
                    plot([i-offset, j-offset],[y,y],'k-');
                    text(mean([i-offset, j-offset]), y*1.05, star_string(p), ...
                        'HorizontalAlignment','center', 'HandleVisibility', 'off');
                end
                % NC stars
                p = pvals_nc{h}(i,j);
                if ~isnan(p) && p < 0.05
                    y = max([mu_nc(i)+sem_nc_vec(i), mu_nc(j)+sem_nc_vec(j)]) * 1.2;
                    plot([i+offset, j+offset],[y,y],'k-');
                    text(mean([i+offset, j+offset]), y*1.05, star_string(p), ...
                        'HorizontalAlignment','center', 'HandleVisibility', 'off');
                end
            end
        end
    end
end


function plot_bar_perDyad_withSEM(d_avg_slp, d_avg_nc, d_sem_slp, d_sem_nc, task_labels, hem_labels, pilot_labels)
% d_avg_slp/nc: cell arrays {1×nTasks}, each element [nDyads × 2]
% d_sem_slp/nc: same structure, SEM across ROIs per dyad × hemisphere
% task_labels:  e.g. {'JS','IR','JR','PN'}
% hem_labels:   e.g. {'Left','Right'}

    nTasks = numel(task_labels);
    nDyads = size(d_avg_slp{1},1);  
    nHems  = size(d_avg_slp{1},2);

    for h = 1:nHems
        figure('Name',[hem_labels{h} ' Hemisphere']); 

        for d = 1:nDyads % loop thru dyad pairs
            mu_slp  = nan(1,nTasks);
            mu_nc   = nan(1,nTasks);
            sem_slp = nan(1,nTasks);
            sem_nc  = nan(1,nTasks);

            for t = 1:nTasks
                mu_slp(t)  = d_avg_slp{t}(d,h);
                mu_nc(t)   = d_avg_nc{t}(d,h);
                sem_slp(t) = d_sem_slp{t}(d,h);
                sem_nc(t)  = d_sem_nc{t}(d,h);
            end

            subplot(ceil(sqrt(nDyads)), ceil(sqrt(nDyads)), d);

            % Bar plot
            bh = bar([mu_slp; mu_nc]','grouped'); 
            hold on;

            % Error bars
            ngroups = nTasks; nbars = 2;
            groupwidth = min(0.8, nbars/(nbars+1.5));
            for i = 1:nbars
                x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
                if i == 1
                    errorbar(x, mu_slp, sem_slp, 'k', 'linestyle','none');
                else
                    errorbar(x, mu_nc, sem_nc, 'k', 'linestyle','none');
                end
            end

            set(gca,'XTick',1:nTasks,'XTickLabel',task_labels);
            title(['Pilot ' pilot_labels(d)]);
            ylabel('Mean ± SEM across ROIs');
            legend({'PWA-clinician','PWA-stranger'},'Location','bestoutside');
        end

        sgtitle([hem_labels{h} ' Hemisphere']); 
    end
end


%% compare SLP-NC for each task incl left and right hems
function X = extract_task_matrix(IBS_struct, dyad_keys, task_field)
% IBS_struct: e.g. IBS_PWA_SLP or IBS_PWA_NC
% dyad_keys: {'dyadic01','dyadic03',...}
% task_field: e.g. 'JS_01','IR_01','JR_01','PN_01'
% Output: [nDyads × 10] matrix (columns = Region1..Region10)

    nD = numel(dyad_keys);
    nR = 10;
    X  = nan(nD, nR);

    for d = 1:nD
        dk = dyad_keys{d};
        for r = 1:nR
            fieldname = "Region" + r;
            if isfield(IBS_struct.(dk).(task_field), fieldname)
                X(d,r) = IBS_struct.(dk).(task_field).(fieldname);
            end
        end
    end
end

% plot
function plot_SLPvsNC_byROI_hemi(X_slp, X_nc, roi_label, task_label, alpha)
% X_slp, X_nc: [nDyads × 10] (columns Region1..Region10)
% roi_label: e.g. 'TPJ'
% task_label: e.g. 'Joint Singing'
% alpha: threshold for sig star (default 0.05)

    if nargin < 5, alpha = 0.05; end

    % Left vs Right columns for this ROI
    roi_idx = find(strcmpi(roi_label, {'IFG','MFG','DLPFC','TPJ','SMG'}));
    colL = roi_idx;         % Regions 15 = Left
    colR = roi_idx + 5;     % Regions 610 = Right

    mu_slp = [mean(X_slp(:,colL),'omitnan'), mean(X_slp(:,colR),'omitnan')];
    mu_nc  = [mean(X_nc(:,colL), 'omitnan'), mean(X_nc(:,colR), 'omitnan')];

    sem_slp = [std(X_slp(:,colL),'omitnan')/sqrt(nnz(~isnan(X_slp(:,colL)))), ...
               std(X_slp(:,colR),'omitnan')/sqrt(nnz(~isnan(X_slp(:,colR))))];
    sem_nc  = [std(X_nc(:,colL),'omitnan')/sqrt(nnz(~isnan(X_nc(:,colL)))), ...
               std(X_nc(:,colR),'omitnan')/sqrt(nnz(~isnan(X_nc(:,colR))))];

    % SLP vs NC test within each hemi
    pL = paired_or_welch(X_slp(:,colL), X_nc(:,colL));
    pR = paired_or_welch(X_slp(:,colR), X_nc(:,colR));

    % ---- plot ----
    figure('Color','w');
    b = bar([mu_slp; mu_nc]'); hold on
    errorbar(b(1).XEndPoints, mu_slp, sem_slp, 'k.', 'LineWidth', 1.2);
    errorbar(b(2).XEndPoints, mu_nc,  sem_nc,  'k.', 'LineWidth', 1.2);

    set(gca,'XTick',[1 2],'XTickLabel',{'Left','Right'}, 'FontSize', 14);
    ylabel('Interbrain Synchrony', 'FontSize', 14);
    legend({'PWA-clinician','PWA-stranger'}, 'Location','northwest', 'FontSize', 10);
    title(sprintf('%s - %s', roi_label, task_label));
    box off; grid off;
    ylim([0 0.6])

    % stars
    p_arr = [pL pR];
    for h = 1:2
        p = p_arr(h);
        x1 = b(1).XEndPoints(h);
        x2 = b(2).XEndPoints(h);
        y  = max([mu_slp(h)+sem_slp(h), mu_nc(h)+sem_nc(h)]) + 0.02;

        if isfinite(p) && p < alpha
            plot([x1 x2],[y y],'k-','LineWidth',1, 'HandleVisibility','off');  % connector
            if     p < 0.001, s='***';
            elseif p < 0.01,  s='**';
            else,             s='*';
            end
            text(mean([x1 x2]), y+0.01, s, ...
                 'HorizontalAlignment','center','FontSize',12);
        end
    end

    % adjust y-limits
    ytop = max([mu_slp+sem_slp, mu_nc+sem_nc]) + 0.08;
    yl = ylim; ylim([min(0,yl(1)) max(yl(2), ytop)]);
    hold off
end

function p = paired_or_welch(x1,x2)
    v1=isfinite(x1); v2=isfinite(x2);
    a=x1(v1); b=x2(v2);
    if numel(a)==numel(b) && numel(a)>=2
        disp("paired")
        [~,p]=ttest(a,b);   % paired
    elseif nnz(v1)>=2 && nnz(v2)>=2
        [~,p]=ttest2(x1(v1),x2(v2),'Vartype','unequal');
    else
        p=NaN;
    end
end


%%

function stats = rm_task_tests(X, task_labels)
% X: [nDyads × nTasks] dyad-level values for one ROI×hemi
% task_labels: {'JS','IR','JR','PN'}
% Returns:
%   stats.mu, stats.sem (1×T)
%   stats.anova_p  (or NaN if fitrm missing)
%   stats.friedman_p
%   stats.p_pair (T×T raw p for paired t-tests)
%   stats.q_pair (T×T FDR-corrected p across all pairs)

    arguments
        X double
        task_labels {mustBeVector}
    end

    [nD, T] = size(X);
    mu  = mean(X, 1, 'omitnan');
    sd  = std (X, 0, 1, 'omitnan');
    N   = sum(isfinite(X), 1);
    sem = sd ./ sqrt(max(N,1));

    % --- Repeated-measures ANOVA if available
    anova_p = NaN;
    try
        tbl = array2table(X, 'VariableNames', matlab.lang.makeValidName(task_labels));
        rm  = fitrm(tbl, sprintf('%s-%s ~ 1', tbl.Properties.VariableNames{1}, tbl.Properties.VariableNames{end}), ...
                    'WithinDesign', table(categorical(task_labels(:)), 'VariableNames', {'Task'}));
        ran = ranova(rm, 'WithinModel','Task');
        anova_p = ran.pValue(1);  % main effect of Task
    catch
        % fitrm/ranova not available; will rely on Friedman
    end

    % --- Friedman (nonparametric RM)
    fr_p = NaN;
    try
        fr_p = friedman(X, 1, 'off');
    catch
    end

    % --- All pairwise paired t-tests (upper triangle)
    p_pair = NaN(T,T);
    for i = 1:T-1
        for j = i+1:T
            xi = X(:,i); xj = X(:,j);
            m  = isfinite(xi) & isfinite(xj);
            xi = xi(m); xj = xj(m);
            if numel(xi) >= 2
                [~,p] = ttest(xi, xj);  % paired across dyads
                p_pair(i,j) = p;
                p_pair(j,i) = p;
            end
        end
    end

    % --- FDR across all unique pairs
    upper = find(triu(true(T),1));
    pvec  = p_pair(upper);
    qvec  = fdr_bh_vector(pvec);
    q_pair = NaN(T,T); q_pair(upper) = qvec; q_pair = q_pair + q_pair.'; % mirror

    stats.mu = mu; stats.sem = sem;
    stats.anova_p = anova_p; stats.friedman_p = fr_p;
    stats.p_pair = p_pair; stats.q_pair = q_pair;
end

function q = fdr_bh_vector(p)
% BenjaminiHochberg FDR for a vector p (ignores NaNs).
    q = nan(size(p));
    v = isfinite(p);
    pv = p(v);
    [ps, idx] = sort(pv);
    m = numel(ps);
    if m==0, q(:) = NaN; return; end
    qtmp = ps .* m ./ (1:m);
    % enforce monotonicity
    for k = m-1:-1:1
        qtmp(k) = min(qtmp(k), qtmp(k+1));
    end
    % scatter back
    q(v) = qtmp(invert_index(idx, m));
    q(q>1) = 1;
end

function o = invert_index(idx, m)
    o = zeros(size(idx)); o(idx) = 1:m;
end

function plot_tasks_SLP_NC_with_two_group_sigs(stats_slp, stats_nc, task_labels, title_nm, varargin)
% stats_slp / stats_nc: from rm_task_tests (must have .mu, .sem, .p_pair, .q_pair)
% task_labels: e.g. {'JS','IR','JR','PN'}
% Options:
%   'UseFDR'       true|false  (default true)
%   'Alpha'        0.05        (default 0.05)
%   'OnlyAdjacent' true|false  (default false)
%   'MaxPairsEach' integer     (default inf)
%   'Pad'          0.02        baseline pad above bars
%   'Step'         0.015       stacking step between lines
%   'YOffsetNC'    0.015       extra vertical offset for NC layer

p = inputParser;
addParameter(p,'UseFDR',false);
addParameter(p,'Alpha',0.05);
addParameter(p,'OnlyAdjacent',false);
addParameter(p,'MaxPairsEach',inf);
addParameter(p,'Pad',0.02);
addParameter(p,'Step',0.015);
addParameter(p,'YOffsetNC',0.015);
parse(p,varargin{:});
opt = p.Results;

mu1  = stats_slp.mu;  sem1 = stats_slp.sem;
mu2  = stats_nc.mu;   sem2 = stats_nc.sem;

pairM_slp = stats_slp.q_pair; 
pairM_nc  = stats_nc.q_pair;
if ~opt.UseFDR
    pairM_slp = stats_slp.p_pair;
    pairM_nc  = stats_nc.p_pair;
end

T = numel(mu1);
figure('Color','w');
b = bar([mu1; mu2]'); hold on
errorbar(b(1).XEndPoints, mu1, sem1, 'k.', 'LineWidth',1.2);
errorbar(b(2).XEndPoints, mu2, sem2, 'k.', 'LineWidth',1.2);
set(gca,'XTick',1:T,'XTickLabel',task_labels, 'FontSize', 12);
ylabel('Interbrain Synchrony', 'FontSize', 14);
legend({'PWA-clinician','PWA-stranger'}, 'Location','northwest', 'FontSize', 10);
title(title_nm); box off; grid off;
ylim([0 0.5])

% --- centers of each task's bar cluster (midpoint of SLP & NC bars) ---
xCenters = mean([b(1).XEndPoints; b(2).XEndPoints], 1);


% Pairs to consider
pairs = [];
if opt.OnlyAdjacent
    for i=1:T-1, pairs(end+1,:) = [i i+1]; end %#ok<AGROW>
else
    for i=1:T-1, for j=i+1:T, pairs(end+1,:)=[i j]; end, end %#ok<AGROW>
end

% Baseline above both groups
y_base = max([mu1+sem1, mu2+sem2]) + opt.Pad;

% ---- Draw SLP sigs (black) ----
draw_sig_layer(xCenters, pairs, pairM_slp, [0, 0.4470, 0.7410], y_base, opt.Step, opt.Alpha, opt.MaxPairsEach);

% ---- Draw NC sigs (gray), offset upward ----
draw_sig_layer(xCenters, pairs, pairM_nc, [0.8500, 0.3250, 0.0980], y_base+opt.YOffsetNC, opt.Step, opt.Alpha, opt.MaxPairsEach);

% Legend entries for lines
plot(nan,nan,'-', 'Color', [0, 0.4470, 0.7410], 'LineWidth',1);      % SLP connector
plot(nan,nan,'-','Color',[0.8500, 0.3250, 0.0980],'LineWidth',1); % NC connector
legend({'PWA-clinician','PWA-stranger'}, 'Location','northwest', 'FontSize', 10);

hold off
end

function draw_sig_layer(xCenters, pairs, pairM, colorSpec, y_base, ystep, alpha, maxPairs)
% Collect significant pairs
hits = [];
for k=1:size(pairs,1)
    i=pairs(k,1); j=pairs(k,2);
    p = pairM(i,j);
    if isfinite(p) && p < alpha
        hits(end+1,:) = [i j p]; %#ok<AGROW>
    end
end
if isempty(hits), return; end
[~,ord]=sort(hits(:,3)); hits=hits(ord,:);
if isfinite(maxPairs) && size(hits,1)>maxPairs
    hits=hits(1:maxPairs,:);
end

% Stack lines to avoid overlap within this layer
T = max(pairs,[],'all');
crossCount = zeros(1,T);

for h=1:size(hits,1)
    i=hits(h,1); j=hits(h,2); p=hits(h,3);
    span = i:j;
    level = max(crossCount(span)) + 1;
    crossCount(span) = level;
    y = y_base + (level-1)*ystep;

    % connector
%     plot([i j],[y y],'-','Color',colorSpec,'LineWidth',1);
%     plot([i i],[y-0.005 y+0.005],'-','Color',colorSpec,'LineWidth',1);
%     plot([j j],[y-0.005 y+0.005],'-','Color',colorSpec,'LineWidth',1);
%     
%     xCenters = mean([b(1).XEndPoints; b(2).XEndPoints],1);
    x1 = xCenters(i);
    x2 = xCenters(j);
    
    plot([x1 x2],[y y],'-','Color',colorSpec,'LineWidth',1);
    plot([x1 x1],[y-0.005 y+0.005],'-','Color',colorSpec,'LineWidth',1);
    plot([x2 x2],[y-0.005 y+0.005],'-','Color',colorSpec,'LineWidth',1);


    % stars
    if     p < 0.001, s='***';
    elseif p < 0.01,  s='**';
    else,             s='*';
    end
    text(mean([x1 x2]), y+0.007, s, ...
     'HorizontalAlignment','center','FontSize',12,'Color',colorSpec);
    %text((i+j)/2, y+0.007, s, 'HorizontalAlignment','center','FontSize',12,'Color',colorSpec);
end
end


function plot_tasks_for_roi(mu, sem, task_labels, title_nm, pair_p, roi_idx, save_dir, ROI_nms, varargin)
% mu, sem: 1×T vectors (means & SEM for each task)
% task_labels: {'JS','IR','JR','PN'}
% title_nm: figure title
% pair_p:  T×T matrix of pairwise p-values (FDR-corrected preferred).
% Optional name/value:
%   'Alpha'        (default 0.05) threshold for significance
%   'OnlyAdjacent' (default false) if true, only test adjacent bars (1-2,2-3,3-4)
%   'MaxPairs'     (default inf)   cap how many sig pairs to draw (avoid clutter)
%   'Pad'          (default 0.02)  vertical pad for first line
%   'Step'         (default 0.015) vertical step between stacked lines

p = inputParser;
addParameter(p,'Alpha',0.05);
addParameter(p,'OnlyAdjacent',false);
addParameter(p,'MaxPairs',inf);
addParameter(p,'Pad',0.02);
addParameter(p,'Step',0.015);
parse(p,varargin{:});
alpha = p.Results.Alpha;
onlyAdj = p.Results.OnlyAdjacent;
maxPairs = p.Results.MaxPairs;
ypad = p.Results.Pad;
ystep = p.Results.Step;

T = numel(mu);

% ---- Base bar plot ----
f = figure('Color','w');
b = bar(1:T, mu); hold on
errorbar(1:T, mu, sem, 'k.', 'LineWidth', 1.2);
set(gca,'XTick',1:T,'XTickLabel',task_labels, 'FontSize', 12);
ylabel('Interbrain Synchrony', 'FontSize', 14);
title(title_nm, 'FontSize', 14);
box off; grid off;

% --- choose which pairs to draw ---
pairs = [];
if onlyAdj
    for i = 1:T-1
        pairs(end+1,:) = [i i+1]; %#ok<AGROW>
    end
else
    for i = 1:T-1
        for j = i+1:T
            pairs(end+1,:) = [i j]; %#ok<AGROW>
        end
    end
end

% Compute a baseline y above all bars for stacking lines
y_base = max(mu + sem) + ypad;

% Collect significant pairs (by p<alpha), sorted from smallest p to largest
P = [];
for k = 1:size(pairs,1)
    i = pairs(k,1); j = pairs(k,2);
    pval = pair_p(i,j);
    if isfinite(pval) && pval < alpha
        P(end+1,:) = [i, j, pval]; %#ok<AGROW>
    end
end
if isempty(P)
    hold off; return; % nothing to draw
end

% Sort by p ascending (strongest first), and cap at MaxPairs
[~,ord] = sort(P(:,3)); P = P(ord,:);
if isfinite(maxPairs) && size(P,1) > maxPairs
    P = P(1:maxPairs,:);
end

% Track how many lines already drawn that cross over each x-position
crossCount = zeros(1,T);

% Draw connectors + stars
for idx = 1:size(P,1)
    i = P(idx,1); j = P(idx,2); pval = P(idx,3);

    % Find how many existing lines cross between i..j to stack this one higher
    span = i:j;
    hlevel = max(crossCount(span)) + 1;
    crossCount(span) = hlevel;

    y = y_base + (hlevel-1)*ystep;

    % horizontal connector with small vertical ticks
    plot([i j], [y y], 'k-', 'LineWidth', 1);              % main line
    plot([i i], [y-0.005 y+0.005], 'k-', 'LineWidth',1);   % left tick
    plot([j j], [y-0.005 y+0.005], 'k-', 'LineWidth',1);   % right tick

    % stars
    if     pval < 0.001, stars = '***';
    elseif pval < 0.01,  stars = '**';
    else                  stars = '*';
    end
    text((i+j)/2, y + 0.007, stars, 'HorizontalAlignment','center', 'FontSize', 12);
end

hold off

roi_name = ROI_nms{roi_idx};
fil_nm_tmp = roi_name;
if contains(title_nm, "PWA-clinician") % tf will be logical 1"PWA-clinician"
    fil_nm_tmp = fil_nm_tmp + "_SLP_";
    if contains(title_nm, "Left")
        fil_nm_tmp = fil_nm_tmp + "LH";
    elseif contains(title_nm, "Right")
        fil_nm_tmp = fil_nm_tmp + "RH";
    end
elseif contains(title_nm, "PWA-stranger")
    fil_nm_tmp = fil_nm_tmp + "_NC_";
    if contains(title_nm, "Left")
        fil_nm_tmp = fil_nm_tmp + "LH";
    elseif contains(title_nm, "Right")
        fil_nm_tmp = fil_nm_tmp + "RH";
    end
end

% save
fig_dir = fullfile(save_dir, "figs");
png_dir = fullfile(save_dir, "png");
if ~isfolder(fig_dir) % Check if the folder does not exist
    mkdir(fig_dir);   % Create the folder
end
if ~isfolder(png_dir) % Check if the folder does not exist
    mkdir(png_dir);   % Create the folder
end
saveas(f, fullfile(fig_dir, fil_nm_tmp + ".fig"));
saveas(f, fullfile(png_dir,  fil_nm_tmp + ".png"));
     
end
% 
% 
% function plot_tasks_for_roi(mu, sem, task_labels, title_nm, pair_q)
% % mu, sem: 1×T
% % pair_q:  T×T FDR-corrected p-values (NaN off-diagonal if not computed)
% % Draw bars with SEM; print significant pairs in Command Window.
% 
%     T = numel(mu);
%     figure('Color','w'); 
%     b = bar(1:T, mu); hold on
%     errorbar(1:T, mu, sem, 'k.', 'LineWidth', 1.2);
%     set(gca,'XTick',1:T,'XTickLabel',task_labels);
%     ylabel('Interbrain Synchrony');
%     title(title_nm);
%     box off; grid off; hold off
% 
%     % Print significant task pairs (q<0.05) to console
%     if nargin >= 5 && ~isempty(pair_q)
%         sigPairs = {};
%         for i = 1:T-1
%             for j = i+1:T
%                 q = pair_q(i,j);
%                 if isfinite(q) && q < 0.05
%                     sigPairs{end+1} = sprintf('%s vs %s: q=%.3g', task_labels{i}, task_labels{j}, q); %#ok<AGROW>
%                 end
%             end
%         end
%         if isempty(sigPairs)
%             fprintf('[%s] No significant task differences after FDR.\n', title_nm);
%         else
%             fprintf('[%s] Significant task pairs (FDR q<0.05):\n', title_nm);
%             fprintf('  - %s\n', sigPairs{:});
%         end
%     end
% end


%
%% Functions
%
% % MAKE THIS A FUNC?
% % % No hemisphere collapse
% % 1) Mean & SEM across dyads for each ROI (no hemi collapse)
% function d_avg_alltasks_rois = avg_dyads_allrois(d_JS_pSLP, , cfg_len)
%     for a = 1:5
%         for m = 1:(cfn_len/2)
%             meanvec(m) = d_JS_pSLP{m}(1,a);
%             rJS_pSLP_avg(1,a) = mean(meanvec, 'omitnan');  % gathering all dyads curr roi ibs value
%             % calc sem across dyads
%             n1 = sum(~isnan(meanvec(m)));
%             rJS_sem(1,a) = std(rJS_avg(1,a), 0, 'omitnan') / sqrt(n1);
%         end
%     end
%     for a = 1:5
%         for m = 1:(cfn_len/2)
%             meanvec(m) = d_JS_pNC{m}(2,a);
%             rJS_pSLP_avg(2,a) = mean(meanvec, 'omitnan'); % mean of each ROI across dyads
%             % calc sem across dyads
%             n1 = sum(~isnan(meanvec(m)));
%             rJS_sem(2,a) = std(rJS_avg(2,a), 0, 'omitnan') / sqrt(n1);
%         end
%     end  % gives you a matrix of average of each dyad for each region for curr task
% end

%% Compare tasks across ROIs (but separate for left and right)
function X = roi_task_matrix(IBS, dyad_keys, task_fields, roi_idx, hemi)
% IBS: struct (e.g., IBS_PWA_SLP or IBS_PWA_NC)
% dyad_keys: {'dyadic01', ...} in the correct patient order
% task_fields: {'JS','IR','JR','PN'} (use your exact task names)
% roi_idx: 1..5  (IFG=1, MFG=2, DLPFC=3, TPJ=4, SMG=5)
% hemi: 'L' or 'R'
% Returns X: [nDyads × nTasks] with per-dyad values

    assert(roi_idx>=1 && roi_idx<=5, 'roi_idx must be 1..5');
    reg = roi_idx + (hemi=='R')*5;  % 1..5 (L) or 6..10 (R)
    nD = numel(dyad_keys); 
    nT = numel(task_fields);
    X  = nan(nD, nT);

    for d = 1:nD
        dk = dyad_keys{d};
        for t = 1:nT
            tf = task_fields{t};
            X(d,t) = IBS.(dk).(tf).("Region" + reg);
        end
    end
end


%% global response - compare slp-nc for each hemisphere
function plot_SLPvsNC_across_tasks_left_or_right(mats_slp, mats_nc, task_labels, hemi_col, title_nm)
% mats_slp / mats_nc: cell arrays of task matrices, each nDyads×2 ([L R])
% task_labels: {'JS','IR','JR','PN'} (same length as mats_*)
% hemi_col: 1 for Left, 2 for Right
% title_nm: title string (e.g., 'Left Hemisphere' or 'Right Hemisphere')
% doFDR: true/false, apply FDR across tasks for this hemisphere

    assert(numel(mats_slp)==numel(mats_nc) && numel(mats_slp)==numel(task_labels), 'Inputs size mismatch.');

    T = numel(mats_slp);
    mu_slp = nan(1,T); sem_slp = nan(1,T); n_slp = nan(1,T);
    mu_nc  = nan(1,T); sem_nc  = nan(1,T); n_nc  = nan(1,T);
    pvals  = nan(1,T);

    for t = 1:T
        x_slp = mats_slp{t}(:, hemi_col);
        x_nc  = mats_nc{t}(:, hemi_col);

        % group stats
        [mu_slp(t), sem_slp(t), n_slp(t)] = group_mean_sem(x_slp);
        [mu_nc(t),  sem_nc(t),  n_nc(t)]  = group_mean_sem(x_nc);

        % paired if possible (same length & index-matched), else Welch t-test
        [~, p] = ttest(x_slp, x_nc);
        %[~, pval(t)] = ttest2(x_slp, x_nc, 'Vartype','unequal');
        [pvals(t)] = p;
    end

    p_for_stars = pvals;

    % ==== PLOT ====
    figure('Color','w');
    b = bar([mu_slp; mu_nc]'); hold on
    % errorbars
    errorbar(b(1).XEndPoints, mu_slp, sem_slp, 'k.', 'LineWidth', 1.2);
    errorbar(b(2).XEndPoints, mu_nc,  sem_nc,  'k.', 'LineWidth', 1.2);

    % stars per task
    for t = 1:T
        s = '';
        p = pvals(t);
        x1 = b(1).XEndPoints(t);
        x2 = b(2).XEndPoints(t);
        y  = max([mu_slp(t)+sem_slp(t), mu_nc(t)+sem_nc(t)]) + 0.02;
        % plot significance if any
        if isfinite(p) && p<=0.05
            plot([x1 x2],[y y],'k-','LineWidth',1)
            s = ''; 
        
            if p < 0.001
                s='***';
            elseif p < 0.01 
                s='**';
            else
                s='*';
            end
        end

        if ~isempty(s)
            text(mean([x1 x2]), y+0.01, s, 'HorizontalAlignment','center','FontSize',12);
        end
    end

    set(gca,'XTick',1:T,'XTickLabel',task_labels, 'FontSize', 12);
    ylim([0 0.45])
    ylabel('Interbrain Synchrony', 'FontSize', 14);
    legend({'PWA-clinician','PWA-stranger'}, 'Location','northeast', 'FontSize', 10); 
    title(title_nm, 'FontSize', 14);
    box off; grid off; hold off
end


%% global - compare hemispheres by averaging across rois 
% stats for comparing hemispheres 
function stats = paired_hemi_stats(X)
% Extract & mask valid pairs
L = X(:,1); R = X(:,2);
m = isfinite(L) & isfinite(R);
L = L(m); R = R(m);

N = numel(L);
stats.N = N;

if N < 2
    stats.meanL   = mean(L,'omitnan');
    stats.meanR   = mean(R,'omitnan');
    stats.meanDiff= mean(L-R,'omitnan');
    stats.t       = NaN; stats.df = NaN; stats.p = NaN;
    stats.ci      = [NaN NaN]; stats.dz = NaN;
    return
end

% Paired t-test (Left vs Right)
[~, pval, ci, st] = ttest(L, R, 'Alpha', 0.05, 'Tail', 'both');
stats.t   = st.tstat;
stats.df  = st.df;
stats.p   = pval;
stats.ci  = ci;      % CI for the mean difference (L-R)

end

function plot_hemi_group(X, title_nm, save_dir, filname)
% X: nDyads × 2 [Left Right]
    [mu, sem, N] = group_mean_sem(X);
    stats = paired_hemi_stats(X);   % your earlier function
    
    f = figure('Color','w'); 
    b = bar([1 2], mu); hold on
    errorbar([1 2], mu, sem, 'k.', 'LineWidth', 1.2)

    % Significance line + stars
    y = max(mu+sem) + 0.02;
    if ~isnan(stats.p) && stats.p < 0.05
        plot([1 2], [y y], 'k-', 'LineWidth', 1);
        if     stats.p < 0.001, stars='***';
        elseif stats.p < 0.01,  stars='**';
        else,                   stars='*';
        end
        text(1.5, y+0.01, stars, 'HorizontalAlignment','center','FontSize',12);
    end

    xlim([0.5 2.5]); set(gca,'XTick',[1 2],'XTickLabel',{'Left','Right'})
    ylim([0 0.45]);
    ylabel('Interbrain Synchrony');
    title(title_nm);
    box off; grid off; hold off

    % save
    fig_dir = fullfile(save_dir, "figs");
    png_dir = fullfile(save_dir, "png");
    if ~isfolder(fig_dir) % Check if the folder does not exist
        mkdir(fig_dir);   % Create the folder
    end
    if ~isfolder(png_dir) % Check if the folder does not exist
        mkdir(png_dir);   % Create the folder
    end
    saveas(f, fullfile(fig_dir, filname + ".fig"));
    saveas(f, fullfile(png_dir,  filname + ".png"));
     
end


%% Left vs right funcs
function plot_leftvsright(mu_2x5, sem_2x5, cats, title_nm, pvals_1x5, save_dir, filname) %, roi_labels)
% mu_2x5, sem_2x5: [2 x 5], rows = [Left; Right], cols = ROIs in order
% cats: categorical(["Left-Brain"; "Right-Brain"])
% pvals_1x5: left-vs-right p-values per ROI (1x5)
% roi_labels: {'IFG','MFG','DLPFC','TPJ','SMG'}
    roi_labels = ['IFG', 'MFG', 'DLPFC', 'TPJ', 'SMG'];
    f = figure(); 
    b = bar(cats, mu_2x5);                 % 5 bar series (one per ROI)
    hold on

    % Add SEM error bars per ROI series
    for k = 1:size(mu_2x5,2)               % k = ROI index
        x = b(k).XEndPoints;               % [x_left, x_right] for this ROI
        m = mu_2x5(:,k);                   % [mu_left; mu_right]
        s = sem_2x5(:,k);                  % [sem_left; sem_right]
        errorbar(x, m, s, 'k', 'linestyle','none', 'LineWidth', 1.2);
    end

    % Significance stars between Left and Right for each ROI
    ypad  = 0.015;                         % pad above the taller bar
    ystep = 0.010;                         % extra pad for the stars
    for k = 1:size(mu_2x5,2)
        p = pvals_1x5(k);
        if ~isnan(p) && p < 0.05
            x = b(k).XEndPoints;           % two x positions (L,R) for this ROI
            m = mu_2x5(:,k);
            s = sem_2x5(:,k);

            y = max(m + s) + ypad;         % line height just above errorbars
            plot([x(1) x(2)], [y y], 'k-', 'LineWidth', 1);

            if     p < 0.001, stars = '***';
            elseif p < 0.01,  stars = '**';
            else               stars = '*';
            end
            text(mean(x), y + ystep, stars, 'HorizontalAlignment','center', 'FontSize', 12);
        end
    end

    % Cosmetics
    %legend(roi_labels, 'Location','eastoutside', 'FontSize',10, 'FontWeight','bold');
    legend('IFG', 'MFG', 'DLPFC', 'TPJ', 'SMG', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Interbrain Synchrony', 'FontSize', 14, 'FontWeight', 'bold');
    ylim([0, max(0.6, max(mu_2x5(:)+sem_2x5(:)) + 0.05)]);
    ax = gca;
    ax.XTickLabel = {'Left-Brain','Right-Brain'};
    ax.XAxis.FontSize = 13; ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontSize = 12; ax.YAxis.FontWeight = 'bold';
    title(title_nm, 'FontSize', 15);
    box off; grid off; hold off
    
    fig_dir = fullfile(save_dir, "figs");
    png_dir = fullfile(save_dir, "png");
    if ~isfolder(fig_dir) % Check if the folder does not exist
        mkdir(fig_dir);   % Create the folder
    end
    if ~isfolder(png_dir) % Check if the folder does not exist
        mkdir(png_dir);   % Create the folder
    end
    saveas(f, fullfile(fig_dir, filname + ".fig"));
    saveas(f, fullfile(png_dir,  filname + ".png"));
end

% function plot_leftvsright(avg_task, cats, title_nm)
%     figure()
%     bJS = bar(cats, avg_task);
%     legend('IFG', 'MFG', 'DLPFC', 'TPJ', 'SMG', 'FontSize', 10, 'FontWeight', 'bold');
%     ylabel('Interbrain Synchrony', 'FontSize', 17, 'FontWeight', 'bold')
%     ylim([-0.025,0.6])
%     
%     % Adjust the font size and weight of the x-axis tick labels
%     ax_JS = gca; % Get current axes
%     ax_JS.XTickLabel = {'Left-Brain', 'Right-Brain'}; % Set custom labels
%     ax_JS.XAxis.FontSize = 15; % Adjust font size
%     ax_JS.XAxis.FontWeight = 'bold'; % Make the labels bold
%     
%     ax_JS.YAxis.FontSize = 14; % Adjust font size for y-axis labels
%     ax_JS.YAxis.FontWeight = 'bold'; % Make y-axis labels bold
%     title(title_nm, 'FontSize', 16)
% end

% stats for comparing hemispheres
function hemi_stats = compare_hemispheres(IBS, dyad_keys, task_fields)
% IBS: struct (SLP or NC)
% dyad_keys: cell array of dyad IDs (e.g. {'dyadic01','dyadic03',...})
% task_fields: e.g. {'JS','IR','JR','PN'}
%
% Compares left vs right (Region1 vs Region6, ..., Region5 vs Region10)

nT = numel(task_fields);
nR = 5; % 5 ROI pairs
pairN = numel(dyad_keys);

pvals = nan(nT, nR);
tvals = nan(nT, nR);
dfs   = nan(nT, nR);

for ti = 1:nT
    tname = task_fields{ti};

    for ri = 1:nR
        left_vals  = nan(pairN,1);
        right_vals = nan(pairN,1);

        for i = 1:pairN
            dkey = dyad_keys{i};

            left_vals(i)  = IBS.(dkey).(tname).("Region"+ri);
            right_vals(i) = IBS.(dkey).(tname).("Region"+(ri+5));
        end

        % Paired t-test left vs right across dyads
        [~, p, ~, st] = ttest(left_vals, right_vals);

        pvals(ti,ri) = p;
        tvals(ti,ri) = st.tstat;
        dfs(ti,ri)   = st.df;
    end
end

hemi_stats.pvals = pvals;  % task × ROIpair
hemi_stats.tvals = tvals;
hemi_stats.dfs   = dfs;
hemi_stats.tasks = task_fields;
hemi_stats.roi_pairs = {'IFG','MFG','DLPFC','TPJ','SMG'}; % left vs right
end

%% stats for comparing SLP-NC
function stats = run_paired_tests(IBS_PWA_SLP, IBS_PWA_NC, slp_keys, nc_keys, task_fields, roi_fields)
% slp_keys, nc_keys: cell arrays of dyad names, aligned by patient order
% task_fields: e.g. {'JS','IR','JR','PN'}
% roi_fields:  e.g. {'IFG','MFG','DLPFC','TPJ','SMG'}
% GIVES US stats for comparing dyad groups

nT = numel(task_fields);
nR = numel(roi_fields);
pairN = min(numel(slp_keys), numel(nc_keys));

pvals = nan(nT,nR);
tvals = nan(nT,nR);
dfs   = nan(nT,nR);
ns    = nan(nT,nR);

for ti = 1:nT
    for ri = 1:nR
        slp_vec = nan(pairN,1);
        nc_vec  = nan(pairN,1);

        for i = 1:pairN
            dslp = slp_keys{i};
            dnc  = nc_keys{i};
            slp_vec(i) = IBS_PWA_SLP.(dslp).(task_fields{ti}).(roi_fields{ri});
            nc_vec(i)  = IBS_PWA_NC.(dnc).(task_fields{ti}).(roi_fields{ri});
        end

        % Paired t-test across dyads
        [~,p,~,st] = ttest(slp_vec, nc_vec);
        pvals(ti,ri) = p;
        tvals(ti,ri) = st.tstat;
        dfs(ti,ri)   = st.df;
        ns(ti,ri)    = numel(slp_vec);
    end
end

stats.pvals = pvals;
stats.tvals = tvals;
stats.dfs   = dfs;
stats.ns    = ns;
stats.tasks = task_fields;
stats.rois  = roi_fields;
end


%% Organize results 
% Organize IBS results for a single task
function [d_JS_pSLP, d_JS_pNC] = single_task_IBS_results(task_name, IBS_PWA_SLP, IBS_PWA_NC, cfn_len, dyad_nms)
    % Reorganizes IBS data froma  struct intoa  cell array for all dyads
    % for a single task. 
    % inputs: 
        % task_name: string 
        % cfn_len: num dyads
    % outputs
        % d_JS_pSLP: PWA-SLP dyad IBS results
            %  cell array containing IBS data for each dyad for task
            % of interest
            % for each cell (dyad); top row regions 1:5 (LH), bottom row regions 6-10 (RH), each column reps an ROI 
        % d_JS_pNC: same as above but for PWA-control dyads

    % Preallocate cell array
    d_JS_pSLP = cell(1, cfn_len);
    d_JS_pNC = cell(1, cfn_len);
    
    temp_vec_pSLP = [];
    temp_vec_pNC = [];
    % Loop over dyads
    for d = 1:cfn_len
        temp_vec_pSLP = [];
        temp_vec_pNC = [];
        %dyad_field = sprintf('dyadic%02d', d); % Dynamic field name for dyad
        dyad_field = dyad_nms(d);
        % Preallocate vector to store region values
        temp_vec = zeros(1, 10);
        % Loop over regions
        for i = 1:10
            region_field = sprintf('Region%d', i); % Dynamic field name for region
            if contains(dyad_field, fieldnames(IBS_PWA_SLP))
                temp_vec_pSLP(i) = IBS_PWA_SLP.(dyad_field).(task_name).(region_field);
            elseif contains(dyad_field, fieldnames(IBS_PWA_NC))
                temp_vec_pNC(i) = IBS_PWA_NC.(dyad_field).(task_name).(region_field);
            end
        end
        %break
        % Reshape into 2x5 matrix
        if ~isempty(temp_vec_pSLP)
            d_JS_pSLP{d} = reshape(temp_vec_pSLP, [5, 2])'; 
        end
        if ~isempty(temp_vec_pNC)
            d_JS_pNC{d} = reshape(temp_vec_pNC, [5, 2])';  % top row regions 1:5 (LH), bottom row regions 6-10 (RH), each column reps an ROI 
        end
    end
    d_JS_pSLP(cellfun('isempty', d_JS_pSLP)) = [];  % remove empties
    d_JS_pNC(cellfun('isempty', d_JS_pNC)) = [];

end

%% Avg hemispheres
function d_cells_avg_hem = avg_hems(d_cells)
% d_cells: cell array, each cell is 2 x nROI (rows=hemispheres, cols=ROIs)
% perDyad: nDyads x nROI, hemisphere-averaged per dyad

    nDyads = numel(d_cells);
    nROI = size(d_cells{1}, 2);
    d_cells_avg_hem = nan(nDyads, nROI);

    for m = 1:nDyads
        % average across hemispheres (rows)
        d_cells_avg_hem(m,:) = mean(d_cells{m}, 1, 'omitnan'); % rows are dyads, cols are ROIs
    end
end

%% "Global responses" - Avg across ROIs on left and right hemispheres 

function d_cells_avg_roi_per_hem = avg_rois_per_hem(d_cells)
% d_cells: cell array, each cell is 2 x nROI (rows=hemispheres, cols=ROIs)
% perDyad: nDyads x nROI, hemisphere-averaged per dyad
% gives a col vector of 2 

    nDyads = numel(d_cells);
    nROI = size(d_cells{1}, 2);
    d_cells_avg_hem = nan(nDyads, nROI);

    for m = 1:nDyads
        % average across hemispheres (rows)
        d_cells_avg_roi_per_hem(m,:) = mean(d_cells{m}, 2, 'omitnan'); % rows are dyads, cols are ROIs
    end
end

function [d_cells_avg_sem_roi_per_hem, sem_perDyad] = avg_sem_rois_per_hem(d_cells)
% d_cells: cell array, each cell is 2 x nROI (rows=hemispheres, cols=ROIs)
% perDyad: nDyads x nROI, hemisphere-averaged per dyad
% gives a col vector of 2 

    nDyads = numel(d_cells);
    nROI = size(d_cells{1}, 2);
    d_cells_avg_hem = nan(nDyads, nROI);

    for m = 1:nDyads
        % average across hemispheres (rows)
        d_cells_avg_sem_roi_per_hem(m,:) = mean(d_cells{m}, 2, 'omitnan'); % rows are dyads, cols are ROIs
        N = sum(~isnan(d_cells{m}), 2);   % 2 x 1
        S = std(d_cells{m}, 0, 2, 'omitnan');
        sem = S ./ sqrt(N); 
        sem_perDyad(m, :) = sem;
    end
end



%% Calc group average and sem on HEMISPHERE AVGED data
function [mu, sem, N] = group_mean_sem(perDyad)
% perDyad: nDyads x nROI
% mu:  1 x nROI (mean across dyads)
% sem: 1 x nROI (standard error)
% N:   1 x nROI (non-NaN counts)
mu = mean(perDyad, 1, 'omitnan');
N  = sum(~isnan(perDyad), 1);
S  = std(perDyad, 0, 1, 'omitnan');   % unbiased std (N-1)
sem = S ./ sqrt(N);
sem(N < 2) = NaN;                     % undefined SEM with <2 samples
end


%% calc group avg and sem without collapsing hemispheres
function [mu_slp, sem_slp, N_slp] = group_mean_sem_hems_sep(d_cells)
    A_slp = cat(3, d_cells{:});   % size: 2 x 5 x nDyads
    nDyads_slp = size(A_slp,3);
    mu_slp  = mean(A_slp, 3, 'omitnan');        % 2 x 5
    S_slp   = std (A_slp, 0, 3, 'omitnan');     % 2 x 5
    N_slp   = sum(~isnan(A_slp), 3);            % 2 x 5
    sem_slp = S_slp ./ sqrt(N_slp); sem_slp(N_slp<2) = NaN;
    % mu_slp(1,:) = LH means per ROI, mu_slp(2,:) = RH means
    
end



%%
% % Hemisphere-collapsed mean & SEM of all dyads (paired, recommended)
% function [hemi_avg, hemi_sem] = calc_group_mean_per_ROI_collapse_hems(d_task_pSLP, d_task_pNC)
% 
%     %nSLP = length(fieldnames(IBS_PWA_SLP));
%     %nNC = length(fieldnames(IBS_PWA_NC));
%     % %
%     nSLP = numel(d_task_pSLP);
%     nNC  = numel(d_task_pNC);
%     
%         %
%     % Average the two hemispheres within each dyad first, then across dyads.
%     hemi_avg = nan(2,5);  % row1 = pSLP, row2 = pNC
%     hemi_sem = nan(2,5);
%     
%     for a = 1:5
%         % pSLP
%         perDyad_slp = nan(1, nSLP);
%         for m = 1:nSLP
%             perDyad_slp(m) = mean(d_task_pSLP{m}(:,a), 'omitnan');  % avg L/R within dyad
%         end
%         hemi_avg(1,a) = mean(perDyad_slp, 'omitnan');
%         n1 = sum(~isnan(perDyad_slp));
%         if n1 >= 2
%             hemi_sem(1,a) = std(perDyad_slp, 0, 'omitnan') / sqrt(n1);
%         else
%             hemi_sem(1,a) = NaN;
%         end
%     
%         % pNC
%         perDyad_nc = nan(1, nNC);
%         for m = 1:nNC
%             perDyad_nc(m) = mean(d_task_pNC{m}(:,a), 'omitnan');
%         end
%         hemi_avg(2,a) = mean(perDyad_nc, 'omitnan');
%         n2 = sum(~isnan(perDyad_nc));
%         if n2 >= 2
%             hemi_sem(2,a) = std(perDyad_nc, 0, 'omitnan') / sqrt(n2);
%         else
%             hemi_sem(2,a) = NaN;
%         end
%     end
%     
%     % hemi_avg -- row 1 = SLP, row 2 = NC
% 
% end
%     
