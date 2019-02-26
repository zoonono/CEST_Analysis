%H-Thebe 29.1.2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       g l u c o C E S T                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Study specific parameters 
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
defaults = {'0','20','1','380','yes','8','11','11','insulin glc','+','30','H-Thebe','no'};
answer1 = inputdlg(prompt,dlg_title,num_lines,defaults);

%rename input parameters
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

clear answer1 defaults num_lines dlg_title prompt

%% set directory - create flags - set initial regularization factor

%set directory that includes the scan folders
directory = 'D:\DATA_glucoCEST\omg_ins_glc_gal\20190129_135535_CEST_Thebe_1_2';

%% experiment timeline
for i = 1: input.number_of_sets
    experiment_time(i) = i* str2num(input.scan_duration)/60; % in minutes
end

%% Make 'images_and_figures' directory into experiment folder
mkdir('images_and_figures');
directory_images = strcat(directory,'\images_and_figures');

%% Put the Bruker numbering of folders of scans in order 
% Checks if the bruker numbering of scans is continuous or not;
% if not, the USER has to type the exact order; it is usually in the correct order already...
[scan_numbering] = scan_ordering(input.first_baseline, input.number_of_scans, input.regular_order); 

%% Read anatomical scan
% Stores some of information (only what's needed) from the files:
% 'read_acqp', 'read_method', 'read_reco' and 'read_2dseq'.
% Note: you can select one on the functions and press F1 to see info (including output)
[parameters_acqp] = read_acqp(input.anatomical_scan_folder_number, directory);
[parameters_method] = read_method(input.anatomical_scan_folder_number, directory);
[parameters_reco] = read_reco(input.anatomical_scan_folder_number, directory);     
[parameters_2dseq] =read_2dseq(input.anatomical_scan_folder_number,...
                               directory, parameters_reco, parameters_acqp);

% Saves all parameters from the four functions in 'anatomical_scan' structure
[anatomical_scan] = structcat(parameters_acqp, parameters_method, ...
                              parameters_reco, parameters_2dseq);
% Remove some variables since they are stored in 'anatomical_scan' struct
clear parameters_acqp   parameters_method   anatomical_scan_folder_number 
clear parameters_reco   parameters_2dseq             

%% Create ROI mask from anatomical image
% Creates a ROI mask using the anatomical image, asks the user to difine the ROIs
% Note: you can select one on the functions and press F1 to see info (including output)
[mask_ROI,image_ROI] = ...
      roi_selection(anatomical_scan.size_x, anatomical_scan.size_y,... 
                    anatomical_scan.FOV, anatomical_scan.image_2dseq,...
                    anatomical_scan.slice_thickness, 1, directory_images); 
% .hdr and .img files of the ROIs are saved in folder images_and_figures
masks = mask_ROI;
[anatomical_scan] = structcat(anatomical_scan, image_ROI);
clear image_ROI mask_ROI

%% Read all CEST scans (baselines and post-injection scans)

h = waitbar(0,'Reading CEST scans...');
for ii = 1 : input.number_of_scans %loop over scan folders  
   % reminder that in "input_parameters.m": 
   % number_of_scans = number_of_baselines + number_of_post_inj_scans 
   % so if the number of the post inj. scans is ZERO, then the number of scans
   % is just the number of baselines.
   % Note: you can select one on the functions and press F1 to see info (including output)
   waitbar(ii/input.number_of_scans)
   [parameters_acqp] = read_acqp (scan_numbering(ii),directory);
   [parameters_method] = read_method(scan_numbering(ii), directory);
   [parameters_reco] = read_reco(scan_numbering(ii), directory);                  
   [parameters_2dseq] =read_2dseq(scan_numbering(ii), directory,...
                                  parameters_reco, parameters_acqp);
   % Saves all parameters from the four functions in 'cest_scan' structure.
   [cest_scan(ii)] = structcat(parameters_acqp, parameters_method, ...
                               parameters_reco, parameters_2dseq);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% separate individual ppm sets within a scan of many sets repeats
i_time = 1;
flag = 0;
for ii = 1 : input.number_of_scans
  ppm_number = 0;  
  for jj = 1 : cest_scan(ii).number_of_ppm_values
    if cest_scan(ii).ppm_values(jj) == cest_scan(ii).ppm_values(1) && jj ~= 1    
      ppm_number = 1; 
      ppms(i_time,ppm_number) = cest_scan(ii).ppm_values(jj);
      data_in_ppm_sets(i_time).image_2dseq(:,:,ppm_number) = cest_scan(ii).image_2dseq(:,:,jj); 
      data_in_ppm_sets(i_time).image_in_gray_scale(:,:,ppm_number) = cest_scan(ii).image_in_gray_scale(:,:,jj);
      data_in_ppm_sets(i_time).ppm_values(ppm_number) = cest_scan(ii).ppm_values(jj);
      data_in_ppm_sets(i_time).slope_values(ppm_number) = cest_scan(ii).slope_values(jj);
      i_time = i_time + 1;
    else
      ppm_number = ppm_number + 1;
      ppms(i_time,ppm_number) = cest_scan(ii).ppm_values(jj);
      data_in_ppm_sets(i_time).image_2dseq(:,:,ppm_number) = cest_scan(ii).image_2dseq(:,:,jj); 
      data_in_ppm_sets(i_time).image_in_gray_scale(:,:,ppm_number) = cest_scan(ii).image_in_gray_scale(:,:,jj);

      data_in_ppm_sets(i_time).ppm_values(ppm_number) = cest_scan(ii).ppm_values(jj);
      data_in_ppm_sets(i_time).slope_values(ppm_number) = cest_scan(ii).slope_values(jj);
    end
  end  
  i_time = i_time + 1;
end

for i = 1 : i_time-1
  data_in_ppm_sets(i).ppm_values = ppms(i,:);
  data_in_ppm_sets(i).number_of_ppm_values = length (data_in_ppm_sets(1).ppm_values);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close(h)
clear parameters   parameters_acqp   parameters_method %stored in CEST scan
clear parameters_reco parameters_2dseq                 %stored in CEST scan
 
%% Resize ROI masks
% Resize the sum_ROI mask to match the cest scan size.
% If ROIs > 1, the sum_ROI is (obviously) the sum of all the individual ROIs.
% It uses the dimensions of the first cest scan to resize the mask created from 
% the anatomical image (usually 256x256) to the size of the cest scans (usually 64x64).
[masks.sum_roi_resized, masks.flag_roi_resized] = ...
       resize_mask(anatomical_scan.size_x, anatomical_scan.size_y, ...
                   cest_scan(1).size_x, cest_scan(1).size_y, masks.sum_roi);
                                 
%Resize individual ROIs (when there is more than one ROI)
if masks.number_of_rois > 1
    for kk = 1 : masks.number_of_rois
    [masks.rois_resized(:,:,kk), masks.flag_roi_resized] = ...
           resize_mask(anatomical_scan.size_x, anatomical_scan.size_y, ...
                       cest_scan(1).size_x, cest_scan(1).size_y, masks.rois(:,:,kk));
    end    
end
                                 
if masks.flag_roi_resized == 1
   disp(['The ROI mask has different dimensions than the (first) CEST image, ',...
         'so now there is also a resized mask!']);
end

% % number of non zero (nnz) elements of resized ROI mask
% nnz_sum_roi = nnz (masks.sum_roi_resized);

% %% Mask the data_in_ppm_sets scans with roi masks
% 
% for ii = 1 : length(data_in_ppm_sets)
%   for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
%     %temporarely created variable image_2dseq for the multiplication    
%     image_2dseq (:,:) = data_in_ppm_sets(ii).image_2dseq(:,:,jj);
% 
%     %scans masked only with RESIZED sum_ROI mask
%     masked_scan(ii).sum_roi(:,:,jj) = immultiply(image_2dseq, masks.sum_roi_resized); 
%     %scans masked with individual ROIs mask
%     if  masks.number_of_rois > 1
%       for kk = 1 : masks.number_of_rois
%       masked_scan(ii).rois(:,:,jj,kk) = immultiply(image_2dseq, masks.rois_resized(:,:,kk)); 
%       end
%     end
%    
%     %convert masked cest scans images in gray scale
%     masked_scan(ii).sum_roi_in_gray_scale(:,:,jj) = mat2gray(masked_scan(ii).sum_roi(:,:,jj));
%     if  masks.number_of_rois > 1
%       for kk = 1 : masks.number_of_rois
%         masked_scan(ii).rois_in_gray_scale(:,:,jj,kk) = mat2gray(masked_scan(ii).rois(:,:,jj,kk));
%       end
%     end
%   end
% end
% clear image_2dseq

%% mask image with roi

%choose manually with mask to use
mask = masks.sum_roi_resized;
%non zero pixels of this mask, for a correct calculation of the mean
nnz_masked_sets = nnz(mask);

for ii = 1 : length(data_in_ppm_sets)
  for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
    %temporarely created variable image_temp for the multiplication    
    image_temp(:,:) = data_in_ppm_sets(ii).image_in_gray_scale(:,:,jj);
    %scans masked only with RESIZED sum_ROI mask
    masked_sets(ii).image_in_gray_scale(:,:,jj) = immultiply(image_temp, mask); 
    
    raw_spectra_roi(ii,jj) = sum(sum(masked_sets(ii).image_in_gray_scale(:,:,jj)))/nnz_masked_sets;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             S P E E D                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set the output format of the command window to the short engineering format 
%with compact line spacing (NO NEED for this command really...)
format shortEng

file = 'Thebe_2901.txt';
path = 'D:\DATA_SPEED\Data\H-Thebe\29012019_Thebe';
cd(path); %go to the path of the specified file
fileID = fopen(file);
formatSpec = '%s %f %f %f %f %f %f %f %s';
A = textscan(fileID,formatSpec,'HeaderLines',12); %ignores first 12 lines
% Col 1: time (hh:mm:ss:ms)
% Col 2: time elapsed from sequence start (sec)
% Col 3: Ch 1 amplitude
% Col 4: Ch 2 amplitude
% Col 5: Ratio Ch1/Ch2
% Col 6: Ratio Ch2/Ch1
% Col 7: Ch1 SNR
% Col 8: Ch2 SNR
% Col 9: date
time    = char(A{1,1}); %(hh:mm:ss:ms)
time_elapsed = A{1,2}; %sec
YFP          = A{1,3}; %ch1
CFP          = A{1,4}; %ch2
FLIIP        = A{1,5}; %(ch1/ch2)
Laconic      = A{1,6}; %(ch2/ch1)
YFP_SNR      = A{1,7}; %(SNR Ch1)
CFP_SNR      = A{1,8}; %(SNR Ch2)
date_exprmt_column  = A{1,9}; %(YYYY:M(M):D(D))
fclose(fileID);

%format the variables you just read!
%make date in format D(D)-M(M)-YYYY
date_exprmt = split(date_exprmt_column(1),":");
date_of_experiment = char(join(flipud(date_exprmt),"-"));

number_of_measurements = length(FLIIP);

time(:,9)= '.'; %(hh:mm:ss.FFF)
[~,~,~,H,MN,S] = datevec(time(:,:));
time_of_the_day_in_min = H*60 + MN + S/60;
clear H MN S
start_time_of_the_day = time_of_the_day_in_min(1);
%experiments time in minutes (start time: 0 min)
time_in_min = time_of_the_day_in_min - start_time_of_the_day;

if length(FLIIP) < 30 
mean_bl_FLIIP = mean(FLIIP(1:length(FLIIP))); %takes the mean of all
norm_FLIIP = 100*(FLIIP/mean_bl_FLIIP-1);
else
mean_bl_FLIIP = mean(FLIIP(1:30)); %takes the mean of the 5min baseline
norm_FLIIP = 100*(FLIIP/mean_bl_FLIIP-1);
end
 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             P L O T S                                   %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = figure;
set(f, 'PaperOrientation', 'landscape'); %'PaperSize', [21 29.7],

% SPEED; plot FLIIP vs TIME  

%    yyy = norm_FLIIP; 
   yyy = movmean(norm_FLIIP,100);
   
   sp3 = subplot(2,1,1);
   p3 = plot(time_in_min(1:length(norm_FLIIP)),yyy);
   p3.LineWidth = 1.5;
   p3.Color = [73/255 0/255 188/255];
   xlim ([0 100])
   ylim ([min(yyy) max(yyy)+ 0.05*max(yyy)]) 

   line([10,55], [-2.5,-2.5],'Color', 'black','LineWidth', 1.5) %([x,x],[y,y])
   line([60,100], [-2.5,-2.5],'Color', 'red','LineWidth', 1.5) %([x,x],[y,y])


ax = ancestor(sp3, 'axes');
   ax.YAxis.FontSize = 13;
   ax.XAxis.FontSize = 13;

   x3 = xlabel('Time (min)');
   x3.FontSize = 12;

   y3 = ylabel('intensity (%)');
   y3.FontSize = 12;

   t = title('FLIIP ratio');
   t.FontSize = 13;

% GlucoCEST; plot timeline of 1.2ppm normalized to first baseline
   
f2 = subplot(2,1,1);
signal2 = - 100*(raw_spectra_roi(:,3)-raw_spectra_roi(1,3))/raw_spectra_roi(1,3);
plot(experiment_time(:),signal2)
min2 = min(signal2);
max2 = max(signal2)+ 0.05*max(signal2);
   line([10,11], [min2 + abs(0.05*min2),min2 + abs(0.05*min2)],'Color', 'black','LineWidth', 1.5) %([x,x],[y,y])
   line([25,70], [min2 + abs(0.05*min2),min2 + abs(0.05*min2)],'Color', 'black','LineWidth', 1.5) %([x,x],[y,y])
   line([83,130],[min2 + abs(0.05*min2),min2 + abs(0.05*min2)],'Color', 'black','LineWidth', 1.5) %([x,x],[y,y])
   line([10,10],   [min2,max2],'Color', [0.8 0.8 0.8],'LineWidth', 0.5) %([x,x],[y,y])
   line([110,110], [min2,max2],'Color', [0.8 0.8 0.8],'LineWidth', 0.5) %([x,x],[y,y])
   line([120,120], [min2,max2],'Color', [0.8 0.8 0.8],'LineWidth', 0.5) %([x,x],[y,y])
   line([130,130], [min2,max2],'Color', [0.8 0.8 0.8],'LineWidth', 0.5) %([x,x],[y,y])
   line([140,140], [min2,max2],'Color', [0.8 0.8 0.8],'LineWidth', 0.5) %([x,x],[y,y])
   ylim ([min2 max2])
   xlim ([0 140])
%plot(experiment_time(21:220),raw_spectra_roi(21:220,3)/raw_spectra_roi(21,3))
ax2 = ancestor(f2, 'axes');
   ax2.YAxis.FontSize = 13;
   ax2.XAxis.FontSize = 13;

   x2 = xlabel('Time (min)');
   x2.FontSize = 12;

   y2 = ylabel('glucoCEST intensity (%)');
   y2.FontSize = 12;

   t = title('glucoCEST @1.2ppm');
   t.FontSize = 13;   
 
   
   % SPEED; plot FLIIP vs TIME  


   
 %%
   cd('D:\DATA_glucoCEST\omg_ins_glc_gal');
   print('-bestfit', 'Thebe2901_.pdf','-dpdf')
