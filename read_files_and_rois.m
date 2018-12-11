%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Afroditi Eleftheriou                                          %
% Date:     13.10.2018                                                    %
% Version:  4.0                                                           %
%           Scripts "read_*.m" etc made into functions                    %
%           Removed intensity mask                                        %
%           Removed average baseline second approach                      %
%                                                                         %
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       %
%                 %         Read CEST data        %                       %
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       %
%                                                                         %
% Reads data from ParaVision 6.0.1 and creates user defined ROIs          %
% NOTES: BS = baseline subtraction (seen in nomenclature)                 %
%        ST = saturation transfer (seen in nomenclature)                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all %prevents some errors popping up when some variables are defined already :D

input_parameters; %User interface for input parameters
cd(directory)

%% Make 'images_and_figures' directory into experiment folder
mkdir('images_and_figures');
directory_images = strcat(directory,'\images_and_figures');

%% Put the Bruker numbering of folders of scans in order 
% Checks if the bruker numbering of scans is continuous or not;
% if not, the USER has to type the exact order 
% it is usually in the correct order already...
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

%% separate individual ppm sets within a scan of many sets repeats
i_time = 1;
flag = 0;
for ii = 1 : input.number_of_scans
  ppm_number = 0;  
  for jj = 1 : cest_scan(ii).number_of_ppm_values
    if cest_scan(ii).ppm_values(jj) == cest_scan(ii).ppm_values(1) && jj ~= 1    
      ppm_number = 1; 
      timeline(i_time,ppm_number) = cest_scan(ii).ppm_values(jj);
      data_in_ppm_sets(i_time) = cest_scan(ii);
      i_time = i_time + 1;
    else
      ppm_number = ppm_number + 1;
      timeline(i_time,ppm_number) = cest_scan(ii).ppm_values(jj);
      data_in_ppm_sets(i_time) = cest_scan(ii);
    end
  end  
  i_time = i_time + 1;
end

for i = 1 : i_time-1
  data_in_ppm_sets(i).ppm_values = timeline(i,:);
  data_in_ppm_sets(i).number_of_ppm_values = length (data_in_ppm_sets(1).ppm_values);
end

%%
close(h)
clear parameters   parameters_acqp   parameters_method %stored in CEST scan
clear parameters_reco parameters_2dseq                 %stored in CEST scan
 
% %% Sort  IF NECESSARY ppm values and image_2dseq and image_in_gray_scale
% % these are the MRI images as acquired by the sequences
% 
% MUST STILL PUT THE IF NECESSARY commands... Now it just sorts them...
% for ii = 1 : input.number_of_scans
% %     %puts unsorted ppm values and images to new variables, only for testing!
% %     cest_scan(ii).unsorted_ppm_values          = cest_scan(ii).ppm_values;
% %     cest_scan(ii).unsorted_image_2dseq         = cest_scan(ii).image_2dseq;
% %     cest_scan(ii).unsorted_image_in_gray_scale = cest_scan(ii).image_in_gray_scale;
%     %sorts the ppm values and images
%     [cest_scan(ii).ppm_values, order]  = sort(cest_scan(ii).ppm_values);
%      cest_scan(ii).image_2dseq         = cest_scan(ii).image_2dseq(:,:,order); 
%      cest_scan(ii).image_in_gray_scale = cest_scan(ii).image_in_gray_scale(:,:,order); 
% end

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

% number of non zero (nnz) elements of resized ROI mask
nnz_sum_roi = nnz (masks.sum_roi_resized);

%% Mask the (ppm sorted) cest scans with anatomical and roi masks

for ii = 1 : input.number_of_scans
  for jj = 1 : cest_scan(ii).number_of_ppm_values
    %temporarely created variable image_2dseq for the multiplication    
    image_2dseq (:,:) = cest_scan(ii).image_2dseq(:,:,jj);

    %scans masked only with RESIZED sum_ROI mask
    masked_scan(ii).sum_roi(:,:,jj) = immultiply(image_2dseq, masks.sum_roi_resized); 
    %scans masked with individual ROIs mask
    if  masks.number_of_rois > 1
      for kk = 1 : masks.number_of_rois
      masked_scan(ii).rois(:,:,jj,kk) = immultiply(image_2dseq, masks.rois_resized(:,:,kk)); 
      end
    end
   
    %convert masked cest scans images in gray scale
    masked_scan(ii).sum_roi_in_gray_scale(:,:,jj) = mat2gray(masked_scan(ii).sum_roi(:,:,jj));
    if  masks.number_of_rois > 1
      for kk = 1 : masks.number_of_rois
        masked_scan(ii).rois_in_gray_scale(:,:,jj,kk) = mat2gray(masked_scan(ii).rois(:,:,jj,kk));
      end
    end
  end
end
clear image_2dseq