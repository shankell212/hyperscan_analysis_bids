
% Save all results for each dyad in a single .mat file for IBS and block 1 and 2 

clear
clc
close all
%%
base_dir = "/projectnb/nphfnirs/s/datasets/Hyperscanning_patients_2025/";
deriv_dir = base_dir + "derivatives";

filesAndDirs = dir(deriv_dir);
names = string({filesAndDirs.name});
%disp(names)
dyad_lst = names(startsWith(names, "dyadic"));
%%
IBS_values = {};
IBS_block1 = {};
IBS_block2 = {};
for i = 1:length(dyad_lst)
%     if i == 7
%         fprintf('Skipping the 7th dyad bc they dont have a second sesh.\n');    % PATIENTS
%         continue; % Skip the rest of the code for this iteration and go to the next.
%     end

    path = fullfile(deriv_dir, dyad_lst(i));

    IBS_file_nm = "IBS_results_" + dyad_lst(i) +".mat";
    load(path + "/" + IBS_file_nm)
    
    IBS_values.(dyad_lst(i)) = dyad_IBS_results;
    IBS_block1.(dyad_lst(i)) = dyad_IBS_block1;
    IBS_block2.(dyad_lst(i)) = dyad_IBS_block2;
    
end
%%
cd(deriv_dir)
save("IBS-values.mat", "IBS_values")
save("IBS-values-block1.mat", "IBS_block1")
save("IBS-values-block2.mat", "IBS_block2")
