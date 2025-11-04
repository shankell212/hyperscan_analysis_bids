clear all;
close all;
clc
%%
%add homer 
cd '/projectnb/nphfnirs/ns/Homer3/'
setpaths

%% Load in MNI table mapping channels to regions
cd '/projectnb/nphfnirs/s/datasets/Hyperscanning_patients_2025/code/'

% MNI.xlsx contains brain regions. They are assigned based on the probe
% group labels are specified in MNI.xlsx sheet 'group_label'
T = readtable("MNI.xlsx",'Sheet','corrected'); % for group 7-12,14
% T = readtable("MNI.xlsx",'Sheet','pilot'); % for group 1-6,13

%% Set paths 
baseDirectory = '/projectnb/nphfnirs/s/datasets/Hyperscanning_patients_2025/';

dyad_num =  '06' %'08'%'07'; % CHANGE (dyad number)

p1 = "BUMA197" %"BUMA079"  %"BUMA304"; % CHANGES FOR EACH DYAD PAIR - subject IDs
p2 = "BUNC278" %"BUNC283"  %"BUNC279";

ses = char(join(['dyadic', dyad_num], ''));

%PN_in_2_runs = True; % CHANGE!!!! THIS IS FOR IF YOU DID 2 SEPARATE SNIRF FILES FOR PN 


%
p1_dir = char(join([baseDirectory, p1, '/', 'ses-', ses,  '/nirs/'], ''));
p2_dir = char(join([baseDirectory, p2, '/', 'ses-', ses,  '/nirs/'], ''));

% initiate struct to store data for paired subjects and their tasks
path1 = [];
path2 = [];

tasks = ["JS", "IR", "JR", "PN", "PN"] %, "PN"]%, "PN"] %, "PN"];
runNums_lst = ["01", "01", "01", "01", "02"] %, "02"]%, "02"];
%run_num = '01';

PREPROCESS = 1;

%% Preprocess fNIRS data
cd '/projectnb/nphfnirs/s/datasets/Hyperscanning_patients_2025/code'
if PREPROCESS
    for i = 1:length(tasks)
        %
        % 4 task blocks include sing together, read separately, read together,
        % and picture naming
        
        %runNum = '01';
        currTask = tasks(i)
        runNum = runNums_lst(i)
        %
        %these are for when the dyad has one or more runs for a certain task -
        %create loop here to iterate and look for correct run number 
        
        if currTask == "PN"
            task_duration = 300; % second % 5 min time block
        elseif currTask == "IR" | currTask == "JR" 
            task_duration = 120;
        else
            task_duration = 120; % 2 min time block
        end
        
        % FIXME: Need to make sure we are grabbing the correct channel/ source
        % dectector
        % Load SNIRF file and preprocess data
        fprintf('Running preprocessing for task %s subject: %s \n', currTask, p1) % CHANGE SUBJECT NUMBER  
        
        % subject 1
        path1.path2snirf{i} = char(join([p1_dir, p1, '_ses-', ses, '_task-', currTask, '_run-', runNum, '_nirs.snirf'],''));
     
        disp(path1.path2snirf{i})
    
        if ~isempty(dir(path1.path2snirf{i}))
            snirf = SnirfLoad(path1.path2snirf{i}); % Load current SNIRF file
            snirf.data(1).dataTimeSeries = hmrR_PreprocessIntensity_Negative(snirf.data); % get rid of negative values in raw data
            path1.snirf{i} = snirf;
            path1.time{i} = snirf.data.time;
    
            [path1.mlActAuto{i}, path1.data_yavg{i}, path1.data_yavgstd{i}, path1.nTrials{i}, path1.data_ynew{i}, path1.data_yresid{i}, path1.data_ysum2{i}, path1.beta_blks{i}, path1.yR_blks{i}, path1.hmrstats{i}] = conc(path1.snirf{i},task_duration);
        
        else
            disp(['No SNIRF file found for task ', currTask, ' and subject ', p1]);
            snirf = [];
            snirf.data(1).dataTimeSeries = [];
            path1.snirf{i} = [];
            path1.time{i} = [];
        end
    
        % subject 2
        fprintf('Running preprocessing for task %s subject %s \n', currTask, p2) % CHANGE SUBJECT NUMBER  
        path2.path2snirf{i} = char(join([p2_dir, p2, '_ses-', ses, '_task-', currTask, '_run-', runNum, '_nirs.snirf'],''));
        disp(path2.path2snirf{i})
        if ~isempty(dir(path2.path2snirf{i}))
            snirf2 = SnirfLoad(path2.path2snirf{i}); % Load current SNIRF file
            snirf2.data(1).dataTimeSeries = hmrR_PreprocessIntensity_Negative(snirf2.data);
            path2.snirf{i} = snirf2;
            path2.time{i} = snirf2.data.time;
    
            [path2.mlActAuto{i}, path2.data_yavg{i}, path2.data_yavgstd{i}, path2.nTrials{i}, path2.data_ynew{i}, path2.data_yresid{i}, path2.data_ysum2{i}, path2.beta_blks{i}, path2.yR_blks{i}, path2.hmrstats{i}] = conc(path2.snirf{i},task_duration);
        else
            disp(['No SNIRF file found for task ', currTask, ' and subject ', p2]);
            snirf2 = [];
            snirf2.data(1).dataTimeSeries = [];
            path2.snirf{i} = [];
            path2.time{i} = [];
        end
        
    end
    
    %if PN_in_2_runs:
    
    
    % Save
    results = char(join([baseDirectory, 'derivatives/', ses, '/'], ''));
    if ~exist(results, 'dir')
        mkdir(results)
    end
    cd(results);
    save(char(join([p1,"_data_", ses], '')),"path1");
    save(char(join([p2,"_data_", ses], '')), "path2"); 

else  % load in data
    results = char(join([baseDirectory, 'derivatives/', ses, '/'], ''));
    path1 = char(join([results, p1, "_data_", ses]));
    load(path1)
    path2 = char(join([results, p2, "_data_", ses]));
    load(path1)

end

%% Load in processed fNIRS data
results = char(join([baseDirectory, 'derivatives/', ses, '/'], ''));
cd(results)
load(char(join([p1,"_data_", ses], ''))) 
load(char(join([p2,"_data_", ses], ''))) 

%% Perform WTC analysis
% corrected probe BA
% BA = [1,4,6,8,9,39,40,44,45,46]; 
groupRegion = 1:10;

WTC = [];

table = T;

AR1_est = 0;
% for k= 1:length(tasks) % CHANGE to x number of tasks    % for dyad 3 JS & IR did not complete for 1 group region -- need to do diff regressive model ??
for k= 1
    for i=1:length(groupRegion)
        fprintf('\n calculating WTC for task %s, groupRegion %d \n', tasks(k), groupRegion(i))
% 
        if i == 2 
            AR1_est = 1;
        elseif i == 4
            AR1_est = 1;
        elseif i == 10
            AR1_est = 1;
        else
            AR1_est = 0;
        end
        %  get the average conc time series across channels specified by BA regions and
        %  brain side
        % get the average conc time series across channels specified by brain IBS
        % regions, incl. IFG,MFG,DLPFC,TPJ,supram
        table.Group = table.Var9; % this is the group region column

        % Find cols w/ 0 and replace w/ NaN (for pruned channels)
        p1_hbo = path1.data_ynew{k}.dataTimeSeries(:,1,:);
        cols_with_zeros = find(any(p1_hbo == 0));
        path1.data_ynew{k}.dataTimeSeries(:,1,cols_with_zeros) = NaN;

        p2_hbo = path2.data_ynew{k}.dataTimeSeries(:,1,:);
        cols_with_zeros = find(any(p2_hbo == 0));
        path2.data_ynew{k}.dataTimeSeries(:,1,cols_with_zeros) = NaN;
       

        % put time series and avg in u
        u_subj1_HbO = [path1.data_ynew{k}.time mean(path1.data_ynew{k}.dataTimeSeries(:,1,find(table.Group==groupRegion(i))),3, 'omitnan')]; % [time, mean time series of channels in current region]
        u_subj2_HbO = [path2.data_ynew{k}.time mean(path2.data_ynew{k}.dataTimeSeries(:,1,find(table.Group==groupRegion(i))),3, 'omitnan')];
        
        %
        % calculate coherence using the HbO data
        % check if subj1 and subj2 have nonzero HbO signals
        if mean(u_subj1_HbO(:,2)==zeros(length(u_subj1_HbO),1)) || mean(u_subj2_HbO(:,2)==zeros(length(u_subj2_HbO),1))
            fprintf('Region # %d contains all zeros and therefore further processing is disabled\n', groupRegion(i))
         
         %assert ( mean(u_subj1_HbO(:,2)~=zeros(length(u_subj1_HbO),1)) || mean(u_subj2_HbO(:,2)~=zeros(length(u_subj2_HbO),1)),  fprintf('Region # %d contains all zeros and therefore further processing is disabled\n', groupRegion(i)))
        
        % Check if all channels are pruned (all NaNs) for current region for each
        % subject -- if yes, skip this region
        elseif all(isnan(u_subj1_HbO(:,2)), 'all') || all(isnan(u_subj2_HbO(:,2)), 'all')
            fprintf('Region # %d contains all NaNs and therefore further processing is disabled\n', groupRegion(i))

        else
            WTC.HbO{k,i}.subj1 = u_subj1_HbO; % save the HbO data
            WTC.HbO{k,i}.subj2 = u_subj2_HbO;
            
            % Attempt to handle potential NaN issues
            if any(isnan(u_subj1_HbO(:))) || any(isnan(u_subj2_HbO(:)))
                warning('NaN detected in HbO signals. This may cause errors in WTC.');
            end
            
            if AR1_est == 1
                u1 = u_subj1_HbO(:); u1 = u1 - mean(u1,'omitnan'); u1 = fillmissing(u1,'linear','EndValues','nearest');
                u2 = u_subj2_HbO(:); u2 = u2 - mean(u2,'omitnan'); u2 = fillmissing(u2,'linear','EndValues','nearest');
                
                A1 = arburg(u1,1); rho1 = -A1(2);
                A2 = arburg(u2,1); rho2 = -A2(2);
                AR1 = [min(max(rho1,0.01),0.99), min(max(rho2,0.01),0.99)];

                % calculate average coherence and periods
                %wtc(u_subj1_HbO,u_subj2_HbO)
                [coherence,period] = wtc(u_subj1_HbO,u_subj2_HbO, 'AR1', AR1);

            else
                % calculate average coherence and periods
                %wtc(u_subj1_HbO,u_subj2_HbO)
                [coherence,period] = wtc(u_subj1_HbO,u_subj2_HbO);
        WTC.coherence{k,i} = coherence; % save coherence (Rsq) score
        WTC.period{k,i} = period; % save period
            end
        end

    end
    %Save results 
    results = char(join([baseDirectory, 'derivatives/', ses, '/'], ''));
    
    cd(results);
    save(char(join(["WTC_GroupRegion_", tasks(k), '_run_', runNums_lst(k), ses, ".mat"],'')), "WTC");
end

%%
%Save results 
results = char(join([baseDirectory, 'derivatives/', ses, '/'], ''));

cd(results);
save(char(join(["WTC_GroupRegion_", ses, ".mat"],'')), "WTC");

%% plot WTC for each region during one task

% plot_results = '/path/to/plot_results';
% mkdir(plot_results)
% cd(plot_results);
plot_dir = char(join([baseDirectory, 'derivatives/', 'plots/', ses, '/'], ''));
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir)
end
cd(plot_dir);

task = {'joint singing', 'separate reading', 'joint reading','picture naming', 'convo'};

% go through all WTC plot
for k=1:length(task) % 4 tasks
    for i=1:length(groupRegion) % 6:length(BA)
        % check if the array is empty
        figure('color',[1 1 1])
        if isempty(WTC.HbO{k,i})
            i = i+1;
        else
            wtc(WTC.HbO{k,i}.subj1,WTC.HbO{k,i}.subj2) % why r we doing this again?
            colorbar('Ticks',[0,0.5,1],'FontSize',24,'FontWeight','bold')
            label = ['Region ' num2str(i) ' WTC during ' task{k}];
            % label = ['R Wernicke Area WTC during ' task{k}];
            title(label,'FontSize',20 )
            ylabel('Period (s)', 'FontSize',20)
            xlabel('Time (s)', 'FontSize',20)

            hold on

            time_arr = [path1.snirf{1,k}.stim(1,1).data(2:end,1) path1.snirf{1,k}.stim(1,2).data(:,1)];
            xline(sort(reshape(time_arr,[1,4])), 'Color',[1 1 1],'LineWidth',5, 'LineStyle','--')
            time_total = path1.time{1,k}(end);

            time_pos = sort(reshape(time_arr,[1,4]))/time_total;


            hold off
            saveas(gcf,[label '.png']);
        end
    end
end


%% Preprocessing functions 
function [mlActAuto, data_yavg, data_yavgstd, nTrials, data_ynew, data_yresid, data_ysum2, beta_blks, yR_blks, hmrstats] = conc(snirf,task_duration)
    %% look at prune channels
    % fs = 1/(wzh.data.time(2) - wzh.data.time(1));
    mlActAuto = hmrR_PruneChannels(snirf.data, snirf.probe, [], [], [0, 1e7], 3, [0.0, 45.0]);
    
    mlActive = mlActAuto{1};
    % wrong probe

    % corrected probe
    index = mlActive(:,2)==[16,17,18,19]; % pilot study
    % index = mlActive(:,2)==[16,18,19,20]; % find the source-detector pairs need to be manually eliminated
    for i = 1:4 % 4 additional ss channels
        mlActive(index(:,i),3)=0;
    end

    mlActAuto = {};
    mlActAuto{1} = mlActive;

    %% intensity to OD
    % Convert intensity to optical density
    dod = hmrR_Intensity2OD(snirf.data);
    %% motion correction
    %dod = hmrR_MotionCorrectSplineSG(dod, mlActAuto, 0.99, 10, 1);
    dod = hmrR_MotionCorrectWavelet(dod, [], mlActAuto, 1.5, 1);
    %% Bandpass filter
    dodBPFilt = hmrR_BandpassFilt(dod, 0.01, 0.5); 
    %0.5 cardiac
    %% OD to concentraion
    % Convert optical density to concentration
    data_y_conc = hmrR_OD2Conc(dodBPFilt, snirf.probe, [1 1]);
   %% GLM
    % Run GLM analysis on the data
    [data_yavg, data_yavgstd, nTrials, data_ynew, data_yresid, data_ysum2, beta_blks, yR_blks, hmrstats] = ...
    hmrR_GLM_CORRECT(data_y_conc, snirf.stim(2), snirf.probe, mlActAuto, [], [], [], [-2.0, task_duration], 1, 1, [1.0 1.0], 15.0, 1, 1, 0); % 15, 2, 3, 0 ?
    
    %hmrR_GLM(data_y_conc, snirf.stim(2), snirf.probe, mlActAuto, [], [], [], [-2.0, task_duration], 1, 1, [1.0 1.0], 15.0, 1, 1, 0); % 15, 2, 3, 0 ?
    %hmrR_GLM(data_y=data_y_conc, stim=snirf.stim(2), probe=snirf.probe, mlActAuto=mlActAuto, Aaux=[], tIncAuto=[], rcMap=[], trange=[-2.0, task_duration+5], glmSolveMethod=1, idxBasis=1, paramsBasis=[1.0 1.0], rhoSD_ssThresh=15.0, flagNuisanceRMethod=1, driftOrder=1, c_vector=0)

end
