%% Batch Analysis of glucoCEST and SPEED data 
% LOG
% created by:   Afroditi     26.02.2019

%add CEST_Analysis to path
addpath('C:\Users\Afroditi\Documents\GitHub\CEST_Analysis');
addpath(genpath('C:\Users\Afroditi\Documents\GitHub\CEST_Analysis'));

%% Load the excel files < with the experiments' information
fname_CEST  = 'Scoresheet_CEST.xlsx';
% fname_SPEED = 'Scoresheet_SPEED.xlsx';

%get CEST and SPEED parameters from excel scoresheet files:
dir_scoresheet_CEST =  'D:\DATA_glucoCEST\Data_organization'; %data from laptop
% dir_scoresheet_SPEED = 'D:\DATA_SPEED'; %data from laptop
ImPath = 1; % 0 = Give image path manually; 1 = Read excel file

%% Read image path 
if  ImPath == 1 %read image path from excel file
    %read glucoCEST excel
    cd(dir_scoresheet_CEST) 
    CEST = readtable(fname_CEST);
    n_cest = size(CEST,1); %number of experiments
       
elseif  ImPath == 0 % Give image path manually
    %eventually write code for setting image path manually 
end

if  n_cest == 0
    error('no experiments in CEST dataset')
end

%% reading the files
for i = 1:n_cest
clearvars -except i n_cest n_speed CEST SPEED;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 load  g l u c o C E S T                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
if CEST{i,15} == 0 %&& CEST{i,14} == 1   
  %CEST{i,14} == 1 -> 1 = ROI file exists
  %CEST{i,15} == 0 -> 0 = analyze, 1 = already analyzed

  %Study specific parameters
  input.date                          = char(CEST{i,2});
  input.mouse_name                    = char(CEST{i,3});
  input.directory_date                = char(CEST{i,4});
  input.number_of_baselines           = CEST{i,5};
  input.number_of_post_inj_scans      = CEST{i,6};
  input.number_of_baseline_sets       = CEST{i,7};
  input.number_of_post_inj_sets       = CEST{i,8};
  input.regular_order                 = CEST{i,9}; %'yes';
  input.anatomical_scan_folder_number = CEST{i,10};
  input.first_baseline                = CEST{i,11};
  input.first_post_inj_scan           = CEST{i,12};
  input.scan_duration                 = CEST{i,13};
  input.administration                = char(CEST{i,16});
  input.fasted                        = CEST{i,44}; %in hours
  
  scan_duration_min_sec = [floor(input.scan_duration/60) mod(input.scan_duration, 60)];
  
  input.number_of_scans = input.number_of_baselines + ...
      input.number_of_post_inj_scans ;
  input.number_of_sets  = input.number_of_baseline_sets + ...
      input.number_of_post_inj_sets ;
  
  %set directory that includes the scan folders
  animal    = input.mouse_name;
  folder      = input.directory_date;
  dir_CEST = fullfile('D:\DATA_glucoCEST\MRI data\CEST',animal, folder); %data from laptop
  
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      load  S  P  E  E  D                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set the output format of the command window to the short engineering format
%with compact line spacing (NO NEED for this command really...)
  format shortEng
if CEST{i,21} == 1 % checks if SPEED was done together with glucoCEST or not
  file_speed = char(CEST{i,23});
  animal     = input.mouse_name;
  folder     = char(CEST{i,24});
  dir_SPEED  = fullfile('D:\DATA_SPEED',animal, folder); %data from laptop
  
  cd(dir_SPEED); %go to the path of the specified file
  fileID = fopen(file_speed);
  formatSpec = '%s %f %f %f %f %f %f %f %s';
  A = textscan(fileID,formatSpec,'HeaderLines',12); %ignores first 12 lines
  time           = char(A{1,1}); %(hh:mm:ss:ms)
  time_elapsed        = A{1,2}; % time elapsed from sequence start (sec)
  YFP                 = A{1,3}; %ch1
  CFP                 = A{1,4}; %ch2
  FLIIP               = A{1,5}; %(ch1/ch2)
  Laconic             = A{1,6}; %(ch2/ch1)
  YFP_SNR             = A{1,7}; %(SNR Ch1)
  CFP_SNR             = A{1,8}; %(SNR Ch2)
  date_exprmt_column  = A{1,9}; %(YYYY:M(M):D(D))
  fclose(fileID);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               load  T E M P E R A T U R E                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CEST{i,18} == 1 %if a temperature log file exists
  cd([dir_CEST,'\temperature_log']);
  file_temperature = char(CEST{i,19});
  fileID = fopen(file_temperature);
  formatSpec = '%f %f';
  AA = textscan(fileID,formatSpec);
  time_temperature = AA{1,1}*10/60; % in minutes (the log's unit is 10sec...)
  temperature      = AA{1,2}; %temperature

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run CEST_SPEED_analysis.m script                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  CEST_SPEED_analysis  
 
  %%
  cd('D:\DATA_glucoCEST\analysis\batch_analysis');
  unique = datestr(datetime,'ddmmyyyyHHMMSS');
  print('-bestfit', [input.mouse_name,'_',date_of_experiment,'_',unique,'.pdf'],'-dpdf')
   
  end
  end
end