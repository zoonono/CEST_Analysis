%47
%% Create intensity-threshold mask from anatomical image
% Creates a mask using the anatomical image, based on pixel intensity 
% Note: you can select one on the functions and press F1 to see info (including output)
[masks, masked_image] = mask_from_image(anatomical_scan.size_x,...  
                                        anatomical_scan.size_y,... 
                                        anatomical_scan.image_2dseq);
%This functions can also create a mask from a CEST image, if the correct input is set ;)
[anatomical_scan] = structcat(anatomical_scan, masked_image);
clear masked_image

%106
%% Resize intensity-threshold (or 'anatomical') mask
% Resize the mask to match the cest scan size.
% It uses the dimensions of the first cest scan to resize the mask created from 
% the anatomical image (usually 256x256) to the size of the cest scans (usually 64x64).
[masks.intensity_resized, masks.flag_intensity_resized] = ...
       resize_mask(anatomical_scan.size_x, anatomical_scan.size_y, ...
                   cest_scan(1).size_x, cest_scan(1).size_y, masks.intensity);      
if masks.flag_intensity_resized == 1
   disp(['The mask of the anatomical image has different dimensions ',...
         'than the (first) CEST image, so now there is also a resized mask!']);
end
% number of non zero (nnz) elements of resized anatomical mask
nnz_intensity = nnz (masks.intensity_resized);

%145
% number of non zero (nnz) elements of resized sum_ROI + intensity segmented mask
masks.intensity_and_roi_resized = immultiply (masks.sum_roi_resized, masks.intensity_resized);
nnz_intensity_and_roi_resized = nnz (masks.intensity_and_roi_resized);

%149 some lines are already removed
%% Mask the (ppm sorted) cest scans with anatomical and roi masks

for ii = 1 : input.number_of_scans
  for jj = 1 : cest_scan(ii).number_of_ppm_values
    %temporarely created variable image_2dseq for the multiplication    
    image_2dseq (:,:) = cest_scan(ii).image_2dseq(:,:,jj);

    %scans masked only with RESIZED intensity-threshold mask
    masked_scan(ii).intensity(:,:,jj) = immultiply(image_2dseq, masks.intensity_resized); 
    % scans masked with both RESIZED ROI mask and RESIZED intensity-threshold mask
    masked_scan(ii).intensity_and_roi(:,:,jj) = ...
          immultiply(masked_scan(ii).sum_roi(:,:,jj), masks.intensity_resized);
    
    % the word "intensity" here refers to intensity-masked data
    masked_scan(ii).intensity_in_gray_scale(:,:,jj) = mat2gray(masked_scan(ii).intensity(:,:,jj));
    masked_scan(ii).intensity_and_roi_in_gray_scale(:,:,jj) = ...
                      mat2gray(masked_scan(ii).intensity_and_roi(:,:,jj));
  end
end
clear image_2dseq