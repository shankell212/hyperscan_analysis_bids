clear
clc
close all
%% set paths

%%
% Analyze hyperscanning data after WTC

base_dir = "/projectnb/nphfnirs/s/datasets/Hyperscanning_patients_2025/";
deriv_dir = base_dir + "derivatives";

% dyad_num = "05";
% task = "JS";
% run = "01";
% 
% filenm = "WTC_GroupRegion_"+task+"_run_"+run+"dyadic"+dyad_num+".mat";
% 
% load(filenm

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

%% Load results
cd(deriv_dir)
load("IBS-values.mat")
load("IBS-values-block1.mat")
load("IBS-values-block2.mat")

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

% IBS = (block1_avg_coh +block2_avg_coh)/2

%% Pull out PWA-SLP and PWA-NC into their own structs
results = IBS_values;
cfn_len = length(fieldnames(IBS_values)); % number of dyads
dyad_names = fieldnames(results);
IBS_PWA_SLP = struct(); %{};  % patient clinician dyads
IBS_PWA_NC = struct(); %{};  % patient healthy control dyads
for i =1:cfn_len
    curr_dyad = string(dyad_names(i))
    
    if contains(curr_dyad, string(T3.DyadNumber(i)))
        dyad_type = string(T3.DyadType(i));
        if contains(dyad_type, "PWA-SLP")
            IBS_PWA_SLP.(curr_dyad) = IBS_values.(curr_dyad);
        elseif contains(dyad_type, "PWA-NC")
            IBS_PWA_NC.(curr_dyad) = IBS_values.(curr_dyad);  % assign correct dyad data to new struct
        end
    end
end


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

% Joint singing
[d_JS_pSLP, d_JS_pNC] = single_task_IBS_results('JS_01', IBS_PWA_SLP, IBS_PWA_NC, cfn_len);  % reorganize data from struct for curr task
d_avg_hem_JS_pSLP = avg_hems(d_JS_pSLP); % average over the hemispheres
d_avg_hem_JS_pNC = avg_hems(d_JS_pNC); % rows are dyads, cols are ROI
[mean_JS_pSLP, sem_JS_pSLP] = group_mean_sem(d_avg_hem_JS_pSLP);  % calc group avg across dyads
[mean_JS_pNC, sem_JS_pNC] = group_mean_sem(d_avg_hem_JS_pNC);  % calc group avg across dyads
mean_JS = vertcat(mean_JS_pSLP, mean_JS_pNC); % concat into 1 matrix. top row = pSLP, bottom row = pNC
sem_JS = vertcat(sem_JS_pSLP, sem_JS_pNC);


% Ind reading
[d_IR_pSLP, d_IR_pNC] = single_task_IBS_results('IR_01', IBS_PWA_SLP, IBS_PWA_NC, cfn_len);
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
[d_JR_pSLP, d_JR_pNC] = single_task_IBS_results('JR_01', IBS_PWA_SLP, IBS_PWA_NC, cfn_len);
d_avg_hem_JR_pSLP = avg_hems(d_JR_pSLP); % average over the hemispheres
d_avg_hem_JR_pNC = avg_hems(d_JR_pNC); % rows are dyads, cols are ROI
[mean_JR_pSLP, sem_JR_pSLP] = group_mean_sem(d_avg_hem_JR_pSLP);  % calc group avg across dyads
[mean_JR_pNC, sem_JR_pNC] = group_mean_sem(d_avg_hem_JR_pNC);  % calc group avg across dyads
mean_JR = vertcat(mean_JR_pSLP, mean_JR_pNC); % concat into 1 matrix. top row = pSLP, bottom row = pNC
sem_JR = vertcat(sem_JR_pSLP, sem_JR_pNC);
%[mean_JR, sem_JR] = calc_group_mean_per_ROI_collapse_hems(d_JR_pSLP, d_JR_pNC); % calc mean across dyads for each ROI 

% Pic Naming
[d_PN_pSLP, d_PN_pNC] = single_task_IBS_results('PN_01', IBS_PWA_SLP, IBS_PWA_NC, cfn_len);
d_avg_hem_PN_pSLP = avg_hems(d_PN_pSLP); % average over the hemispheres
d_avg_hem_PN_pNC = avg_hems(d_PN_pNC); % rows are dyads, cols are ROI
[mean_PN_pSLP, sem_PN_pSLP] = group_mean_sem(d_avg_hem_PN_pSLP);  % calc group avg across dyads
[mean_PN_pNC, sem_PN_pNC] = group_mean_sem(d_avg_hem_PN_pNC);  % calc group avg across dyads
mean_PN = vertcat(mean_PN_pSLP, mean_PN_pNC); % concat into 1 matrix. top row = pSLP, bottom row = pNC
sem_PN = vertcat(sem_PN_pSLP, sem_PN_pNC);
%[mean_PN, sem_PN] = calc_group_mean_per_ROI_collapse_hems(d_PN_pSLP, d_PN_pNC); % calc mean across dyads for each ROI 

disp('Group average across dyads calculated for each task and ROI')


%% Avg of all dyads for each ROI for each task

mean_task_lst = {mean_JS, mean_IR, mean_JR, mean_PN};
sem_task_lst = {sem_JS, sem_IR, sem_JR, sem_PN};

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

%% Perform Wilcoxon sign ranked test
% Pilot IDs per group, in the same order the group rows appear in T:
slp_ids = T3.pilot( T3.DyadType, contains("SLP"));     % column vector
nc_ids  = T3.pilot( T3.DyadType, contains("NC")); 

% Perform Wilcoxon signed-rank ttest 
% organize dyads by patients


[id_slp,~,g1] = unique(string(id_pSLP));
[id_nc, ~,g2] = unique(string(id_pNC));
[X_aligned, Y_aligned] = deal(nan(numel(intersect(id_slp,id_nc)), nROI));
[common, ia, ib] = intersect(id_slp, id_nc, 'stable');

X_aligned = X(ia,:);  % pSLP per patient
Y_aligned = Y(ib,:);  % pNC per patient














%
%% Functions
%

% MAKE THIS A FUNC?
% % No hemisphere collapse
% 1) Mean & SEM across dyads for each ROI (no hemi collapse)
% for a = 1:5
%     for m = 1:(cfn_len/2)
%         meanvec(m) = d_JS_pSLP{m}(1,a);
%         rJS_pSLP_avg(1,a) = mean(meanvec, 'omitnan');  % gathering all dyads curr roi ibs value
%         % calc sem across dyads
%         n1 = sum(~isnan(meanvec(m)));
%         rJS_sem(1,a) = std(rJS_avg(1,a), 0, 'omitnan') / sqrt(n1);
%     end
% end
% for a = 1:5
%     for m = 1:(cfn_len/2)
%         meanvec(m) = d_JS_pNC{m}(2,a);
%         rJS_pSLP_avg(2,a) = mean(meanvec, 'omitnan'); % mean of each ROI across dyads
%         % calc sem across dyads
%         n1 = sum(~isnan(meanvec(m)));
%         rJS_sem(2,a) = std(rJS_avg(2,a), 0, 'omitnan') / sqrt(n1);
%     end
% end  % gives you a matrix of average of each dyad for each region for curr task
% 


% Organize IBS results for a single task
function [d_JS_pSLP, d_JS_pNC] = single_task_IBS_results(task_name, IBS_PWA_SLP, IBS_PWA_NC, cfn_len)
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
        dyad_field = sprintf('dyadic%02d', d); % Dynamic field name for dyad
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

% Calc group average and sem on HEMISPHERE AVGED data
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


% calc group avg and sem without collapsing hemispheres
function [mu, sem, N] = group_mean_sem_hems_sep(d_cells)
    A_slp = cat(3, d_cells{:});   % size: 2 x 5 x nDyads
    nDyads_slp = size(A_slp,3);
    mu_slp  = mean(A_slp, 3, 'omitnan');        % 2 x 5
    S_slp   = std (A_slp, 0, 3, 'omitnan');     % 2 x 5
    N_slp   = sum(~isnan(A_slp), 3);            % 2 x 5
    sem_slp = S_slp ./ sqrt(N_slp); sem_slp(N_slp<2) = NaN;
    % mu_slp(1,:) = LH means per ROI, mu_slp(2,:) = RH means
    
end


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

