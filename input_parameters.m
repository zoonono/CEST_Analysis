%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Afroditi Eleftheriou                                                    %
% 22.11.2017                                                              %
%                                                                         %
% This part of the code uses a dialog window for the user to put the      %
% initial parameters of the analysis.                                     %
% Also sets the working directory.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Study specific parameters 

prompt = {'Number of baseline scans',...
          'Number of baseline sets',...
          'Number of post-injection scans',...
          'Number of post-injection sets',...
          'Uninterapted numbering of scans (yes/no)',...
          'Anatomical scan number',...
          'First baseline folder number',...
          'First post-injection folder number',...
          'Administration agent/concentration/route',...
          'Positive or negative ppm of agent?',...
          'Actual scan duration (in sec)',...
          'Name of mouse :)',...
          'Fasted or not (yes/no)'};
dlg_title = 'Study specific input';
num_lines = 1;
defaults = {'0','0','1','200','yes','5','10','10','3omg insulin glc gal','+','30','H-Thebe','no'};
% defaults = {'1','20','5','280','yes','5','9','10','3omg insulin glc gal','+','30','H-Thebe','no'};
% defaults = {'1','20','5','280','yes','7','12','13','3omg insulin glc gal','+','30','H-Adrastea2012','no'};
% defaults = {'1','20','5','280','yes','6','9','10','3omg insulin glc gal','+','30','H-Io','no'};
% defaults = {'1','20','5','280','yes','6','9','10','3omg insulin glc gal','+','30','H-Metis','no'};
answer1 = inputdlg(prompt,dlg_title,num_lines,defaults);


%% Analysis initial parameters 

prompt = {'R^2 threshold (''goodness'' of interpolation):',...
          'At what ppm to calculate AUC from STasym?:',...
          'Width of AUC interval:',...
          'For the second AUC calculation, starting from 0 ppm until:',...
          '''custom'' or ''automatic'' regularization factor of interpolation?',...
          'Regularization factor of cubic B-spline interpolation (in vivo: high, e.g. 0.99; in vitro: low, e.g. 0.3):',...
          'Display the maps with ROI mask (type: roi) or morphological mask (type: morph)?'};
dlg_title = 'Initial parameters for analysis';
num_lines = 1;
defaults = {'0.97','1.2','1.0','3.0','custom','0.99','roi'};
answer2 = inputdlg(prompt,dlg_title,num_lines,defaults);

%% rename input parameters

input.number_of_baselines = str2num(answer1{1});
input.number_of_baseline_sets = str2num(answer1{2});
input.number_of_post_inj_scans = str2num(answer1{3});
input.number_of_post_inj_sets = str2num(answer1{4});

input.number_of_scans = input.number_of_baselines + ...
                        input.number_of_post_inj_scans ;

input.number_of_sets = input.number_of_baseline_sets + ...
                        input.number_of_post_inj_sets ;
                    

input.regular_order = answer1{5};
input.anatomical_scan_folder_number = answer1{6};
input.first_baseline = str2num(answer1{7});
input.first_post_inj_scan = str2num(answer1{8});
input.administration = answer1{9};
input.agent_sign = answer1{10};
input.scan_duration = answer1{11};
  scan_duration_min_sec = ...
  [floor(str2num(input.scan_duration)/60) mod(str2double(input.scan_duration), 60)];
input.mouse_name = answer1{12};
input.fasted = answer1{13};

input.Rsqr_threshold = str2num(answer2{1}); % R^2, "goodness" of pixel value threshold!
input.ppm_mid_point_of_auc = str2num(answer2{2});
input.width_integral = str2num(answer2{3});
input.endpoint_auc_from_zero = str2num(answer2{4});
input.auto_or_custom = answer2{5};
input.regularization_weight = str2num(answer2{6});
input.roi_or_anatomical = answer2{7};

clear answer1 answer2 defaults num_lines dlg_title prompt

%% set directory - create flags - set initial regularization factor

%set directory that includes the scan folders
cd ('D:\DATA_glucoCEST');
directory=uigetdir;

%creates flags to be used in the interpolation calculation
if strcmp(input.agent_sign,'+')
input.sign_st = 1;
input.flag_integral = 0;
end
if strcmp(input.agent_sign,'-')
input.sign_st = -1;
input.flag_integral = 1;
end

%% experiment timeline
for i = 1: input.number_of_sets
    experiment_time(i) = i* str2num(input.scan_duration)/60; % in minutes
end
