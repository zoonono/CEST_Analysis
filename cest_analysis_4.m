%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Afroditi Eleftheriou                                          %
% Date:     12.10.2018                                                    %
% Version:  4.0                                                           %
%                                                                         %
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       %
%                 %  Main CEST analysis program!  %                       %
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       %
%                                                                         %
% Analysis of data from ParaVision 6.0.1                                  %
% NOTES: run first "read_files_and_rois.m" or import files                %
%        BS = baseline subtraction (seen in nomenclature)                 %
%        ST = saturation transfer (seen in nomenclature)                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Interpolation of pixel-wise data and R^2 calculation
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%-%         P I X E L - W I S E   I N T E R P O L A T I O N  ! !        %-%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

%interpolate pixel-wise z-spectra (ONLY FOR PIXEL WITHIN sum_roi !!!) 
%Also calculates R^2 of each pixel ("goodness of interpolation")
h = waitbar(0,'Pixel-wise z-spectra interpolation...');
for ii = 1 : input.number_of_scans
  waitbar(ii/input.number_of_scans)
  number_of_pixel = 0;  %within each scan it'll calculate the number of non zero pixels
  for jj = 1 : size(masked_scan(ii).sum_roi,1)
    for kk = 1 : size(masked_scan(ii).sum_roi,2)
      if masks.sum_roi_resized(jj,kk) > 0 %check if the pixel is within the sum_ROI
        number_of_pixel = number_of_pixel + 1;
        %temporary struct name "interpolated", will be named "interpolated_pixels"  
        %the data get NORMALIZED in the "B_spline_Rsqr_STAsym" function
        [interpolated] = B_spline_Rsqr_STAsym(cest_scan(ii).ppm_values, ... 
                         reshape(masked_scan(ii).sum_roi(jj,kk,:), ...
                         [1,cest_scan(ii).number_of_ppm_values]),...
                         input.regularization_weight, input.auto_or_custom,...
                         input.sign_st);
        interpolated.pixel_position = [jj,kk];
        %put all scans in one (renamed) struct, data normilized (done in B_spline_Rsqr_STAsym function)
        interpolated_pixels(ii,number_of_pixel)=interpolated; 
        image_Rsqr(ii,jj,kk) = interpolated_pixels(ii,number_of_pixel).Rsqr_of_z_sprectrum;
        % image_rsqr(number of scan, x, y)
      else
        %if the pixel value is zero, does not need to interpolate anything
        image_Rsqr(ii,jj,kk) = 0.0 ;
      end
    end
  end
end
close(h)
clear interpolated

%% Create Rsqr_threshold MASKs 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for RESIZED sum_roi masked data %
masks.Rsqr_sum_roi_resized = zeros(size(masks.sum_roi_resized));
% Goes through all the non zero pixels and checks the R^2 value of the pixel
% Rejects pixels that don't have a good interpolation curve relative to the raw data, 
%by checking if the pixel's R^2 is greater than or equal to the R^2 threshold value
for ll = 1 : size(interpolated_pixels,2) % 1 : number of non-zero pixels
  if interpolated_pixels(ii,ll).Rsqr_of_z_sprectrum >= input.Rsqr_threshold
    i1 = interpolated_pixels(ii,ll).pixel_position(1); %to avoid a big index name!
    i2 = interpolated_pixels(ii,ll).pixel_position(2); %avoid big index name in the next line!
    masks.Rsqr_sum_roi_resized(i1,i2) = 1;
  end
end
clear i1 i2

%%%%%%%%%%%%%%%%%%%%%%%
% for individual ROIs %
if masks.number_of_rois > 1 
  for kk = 1 : masks.number_of_rois
    % ".*" is for element-wise multiplication (same as immultiply)
    masks.Rsqr_rois_resized(:,:,kk) = masks.Rsqr_sum_roi_resized(:,:) .* masks.rois_resized(:,:,kk);
  end
end

%% Mask data with the Rsqr_threshold mask

for ii = 1 : input.number_of_scans
  for jj = 1 : cest_scan(ii).number_of_ppm_values
    %temporarely created variable image_2dseq for the multiplication    
    image_2dseq (:,:) = cest_scan(ii).image_2dseq(:,:,jj);

    % for RESIZED sum_ROI masked data
    masked_scan(ii).Rsqr_sum_roi(:,:,jj) = immultiply(image_2dseq, masks.Rsqr_sum_roi_resized); 
    % for RESIZED individual ROIs
    if  masks.number_of_rois > 1
      for kk = 1 : masks.number_of_rois
      masked_scan(ii).Rsqr_rois(:,:,jj,kk) = immultiply(image_2dseq, masks.Rsqr_rois_resized(:,:,kk)); 
      end
    end
  
    %convert masked cest scans images in gray scale
    masked_scan(ii).Rsqr_sum_roi(:,:,jj) = mat2gray(masked_scan(ii).Rsqr_sum_roi(:,:,jj));
    if  masks.number_of_rois > 1
      for kk = 1 : masks.number_of_rois
        masked_scan(ii).Rsqr_rois_in_gray_scale(:,:,jj,kk) = mat2gray(masked_scan(ii).Rsqr_rois(:,:,jj,kk));
      end
    end
  end
end
clear image_2dseq

%% mean values and interpolation of MASKED data
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%-%                                                                     %-%
%-%   M E A N   V A L U E S   A  N D   I N T E R P O L A T I O N  ! !   %-%
%-%                                                                     %-%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%%%%%%%%%%%%%%%
% for sum_ROI %
%%%%%%%%%%%%%%%
%mean value per ppm image, over sum_ROI masked pixel values
for ii = 1 : input.number_of_scans
    for jj = 1 : cest_scan(ii).number_of_ppm_values
        z_spectra_raw(ii).mean_sum_roi(jj)= mean2(masked_scan(ii).sum_roi(:,:,jj));
    end
end
%calculate the normalized_mean_sum_roi
for ii = 1 : input.number_of_scans
  max_value_of_roi = max (z_spectra_raw(ii).mean_sum_roi(jj));
  for jj = 1 : cest_scan(ii).number_of_ppm_values
  z_spectra_raw(ii).normalized_sum_roi(jj) = ...
                  z_spectra_raw(ii).mean_sum_roi(jj)/ max_value_of_roi;           
  end
end
%INTERPOLATE z_spectra_raw.mean_sum_roi
for ii = 1: input.number_of_scans
  [interpolated_sum_roi(ii)] = ...  %gets normilized in B_spline_Rsqr_STAsym function
  B_spline_Rsqr_STAsym(cest_scan(ii).ppm_values, z_spectra_raw(ii).mean_sum_roi,...
           input.regularization_weight, input.auto_or_custom,input.sign_st); 
end

%%%%%%%%%%%%%%%%%%%%
% for Rsqr_sum_ROI %
%%%%%%%%%%%%%%%%%%%%
%mean value per ppm image, over Rsqr_sum_ROI masked pixel values
for ii = 1 : input.number_of_scans
    for jj = 1 : cest_scan(ii).number_of_ppm_values
        z_spectra_raw(ii).mean_Rsqr_sum_roi(jj)= mean2(masked_scan(ii).Rsqr_sum_roi(:,:,jj));
    end
end
%calculate the normalized_mean_Rsqr_sum_roi
for ii = 1 : input.number_of_scans
  max_value_of_roi = max (z_spectra_raw(ii).mean_Rsqr_sum_roi(jj));
  for jj = 1 : cest_scan(ii).number_of_ppm_values
  z_spectra_raw(ii).normalized_Rsqr_sum_roi(jj) = ...
                    z_spectra_raw(ii).mean_Rsqr_sum_roi(jj)/ max_value_of_roi;           
  end
end
%INTERPOLATE z_spectra_raw.mean_Rsqr_sum_roi
for ii = 1: input.number_of_scans
  [interpolated_Rsqr_sum_roi(ii)] = ...  %gets normilized in B_spline_Rsqr_STAsym function
  B_spline_Rsqr_STAsym(cest_scan(ii).ppm_values, z_spectra_raw(ii).mean_Rsqr_sum_roi,...
           input.regularization_weight, input.auto_or_custom,input.sign_st); 
end
clear max_value_of_roi

%%%%%%%%%%%%%%%%%%%%%%%
% for individual ROIs %
%%%%%%%%%%%%%%%%%%%%%%%
%mean value per ppm image, over individual ROI masked pixel values
if masks.number_of_rois > 1 
  for ii = 1 : input.number_of_scans
    for kk = 1 : masks.number_of_rois
      for jj = 1 : cest_scan(ii).number_of_ppm_values
        z_spectra_raw(ii).rois(kk,jj)= mean2(masked_scan(ii).rois(:,:,jj,kk)); 
      end
    end
  end
  %calculate the normalized_rois
  for ii = 1 : input.number_of_scans
    for kk = 1 : masks.number_of_rois
      max_values_of_rois(kk) = max (z_spectra_raw(ii).rois(kk,:));
      for jj = 1 : cest_scan(ii).number_of_ppm_values
      z_spectra_raw(ii).normalized_rois(kk,jj) = ...
                        z_spectra_raw(ii).rois(kk,jj)/ max_values_of_rois(kk);
      end
    end
  end
  %INTERPOLATE z_spectra_raw.rois
  for ii = 1: input.number_of_scans
    for kk = 1: masks.number_of_rois
    [interpolated_rois(kk,ii)] = ... %gets normilized in B_spline_Rsqr_STAsym function
    B_spline_Rsqr_STAsym(cest_scan(ii).ppm_values, z_spectra_raw(ii).rois(kk,:),...
             input.regularization_weight, input.auto_or_custom,input.sign_st);
    end
  end 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % for individual Rsqr_ROIs %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%% (already within "if #rois > 1" statement)
  %mean value per ppm image, over individual ROI masked pixel values
  for ii = 1 : input.number_of_scans
    for kk = 1 : masks.number_of_rois
      for jj = 1 : cest_scan(ii).number_of_ppm_values
        z_spectra_raw(ii).Rsqr_rois(kk,jj)= mean2(masked_scan(ii).Rsqr_rois(:,:,jj,kk)); 
      end
    end
  end
  %calculate the normalized_rois
  for ii = 1 : input.number_of_scans
    for kk = 1 : masks.number_of_rois
      max_values_of_rois(kk) = max (z_spectra_raw(ii).Rsqr_rois(kk,:));
      for jj = 1 : cest_scan(ii).number_of_ppm_values
      z_spectra_raw(ii).normalized_Rsqr_rois(kk,jj) = ...
                        z_spectra_raw(ii).Rsqr_rois(kk,jj)/ max_values_of_rois(kk);
      end
    end
  end
  %INTERPOLATE z_spectra_raw.rois
  for ii = 1: input.number_of_scans
    for kk = 1: masks.number_of_rois
    [interpolated_Rsqr_rois(kk,ii)] = ... %gets normilized in B_spline_Rsqr_STAsym function
    B_spline_Rsqr_STAsym(cest_scan(ii).ppm_values, z_spectra_raw(ii).Rsqr_rois(kk,:),...
             input.regularization_weight, input.auto_or_custom,input.sign_st);
    end
  end
end %end of "if masks.number_of_rois > 1" 

%%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%-%                                                                     %-%
%-%          B A S E L I N E   S U B T R A C T I O N   ! ! !            %-%
%-%                                                                     %-%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%  

%% Average baseline z-spectrum
avg_of_raw_baselines_and_interpolation

%% Subtraction of average baseline
baseline_subtraction %(BS)

%% AUC
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%-%                                                                     %-%
%-%          A R E A   U N D E R   T H E   C U R V E   ! ! !            %-%
%-%                                                                     %-%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%done on STAsym interpolated data (interpolated_*.STAsym) and 
%baseline-subtraction interpolated data (BS_interpolated_mean_*)

%% AUC for STAsym interpolated data
for ii = 1 : input.number_of_scans
  % sum_roi
  [auc_STAsym_sum_roi(ii)] = ...
       auc(interpolated_sum_roi(ii).ppm_positive, interpolated_sum_roi(ii).STAsym, ...
           input.ppm_mid_point_of_auc, input.width_integral, ...
           input.endpoint_auc_from_zero);
  % Rsqr_sum_roi
  [auc_STAsym_Rsqr_sum_roi(ii)] = ...
       auc(interpolated_Rsqr_sum_roi(ii).ppm_positive, ...
           interpolated_Rsqr_sum_roi(ii).STAsym, ...
           input.ppm_mid_point_of_auc, input.width_integral, ...
           input.endpoint_auc_from_zero);
  % Rsqr_rois
  if masks.number_of_rois > 1 
    for kk = 1 : masks.number_of_rois
      [auc_STAsym_Rsqr_rois(kk,ii)] = auc(interpolated_Rsqr_rois(kk,ii).ppm_positive, ...
                                        interpolated_Rsqr_rois(kk,ii).STAsym, ...
                                        input.ppm_mid_point_of_auc, input.width_integral, ...
                                        input.endpoint_auc_from_zero);      
    end
  end
  % pixels
  for ll = 1 : size(interpolated_pixels,2)
    [auc_STAsym_pixels(ii,ll)] = auc(interpolated_pixels(ii,ll).ppm_positive, ...
                                     interpolated_pixels(ii,ll).STAsym, ...
                                     input.ppm_mid_point_of_auc, ...
                                     input.width_integral, ...
                                     input.endpoint_auc_from_zero);
  end
end

%% AUC for baseline-subtraction interpolated data
% only done on the 1st approach BS data that is to average the raw baselines 
% and then interpolate ("BS_interpolated_mean_*" variables).
for ii = 1 : input.number_of_scans
  % sum_roi
  [auc_BS_sum_roi(ii)] = auc(interpolated_sum_roi(ii).ppm, ...
                             BS_interpolated_mean_sum_roi(ii,:), ...
                             input.ppm_mid_point_of_auc, input.width_integral, ...
                             input.endpoint_auc_from_zero);
  % Rsqr_sum_roi
  [auc_BS_Rsqr_sum_roi(ii)] = auc(interpolated_Rsqr_sum_roi(ii).ppm, ...
                                BS_interpolated_mean_Rsqr_sum_roi(ii,:), ...
                                input.ppm_mid_point_of_auc, input.width_integral, ...
                                input.endpoint_auc_from_zero);
  % Rsqr_rois
  if masks.number_of_rois > 1 
    for kk = 1 : masks.number_of_rois
      [auc_BS_Rsqr_rois(kk,ii)] = auc(interpolated_Rsqr_rois(kk,ii).ppm, ...
                                    BS_interpolated_mean_Rsqr_rois(kk,ii,:), ...
                                    input.ppm_mid_point_of_auc, input.width_integral, ...
                                    input.endpoint_auc_from_zero);      
    end
  end
  % pixels 
  for ll = 1 : size(interpolated_pixels,2)
    [auc_BS_pixels(ii,ll)] = auc(interpolated_pixels(ii,ll).ppm, ...
                             BS_interpolated_mean_pixels(ii,:,ll), ...
                             input.ppm_mid_point_of_auc, input.width_integral, ...
                             input.endpoint_auc_from_zero);
  end
end

%% Calculating the mean and std of the baselines and post injection scans
%Calculating the mean and std of the baselines 
A = zeros(1,input.number_of_baselines); 
C = zeros(1,input.number_of_baselines); 
for i = 1: input.number_of_baselines
    A(i) = auc_STAsym_sum_roi(i).auc_around_mid_point;
    C(i) = auc_BS_sum_roi(i).auc_around_mid_point;
end
mean_auc_glucose_STAsym_sum_roi = mean(A);
 std_auc_glucose_STAsym_sum_roi = std(A);
mean_auc_glucose_BS_sum_roi     = mean(C);
 std_auc_glucose_BS_sum_roi     = std(C);

%Calculating the mean and std of the post injection scans 
A1 = zeros(1,input.number_of_post_inj_scans); 
C1 = zeros(1,input.number_of_post_inj_scans); 
for i = 1: input.number_of_post_inj_scans
    A1(i) = auc_STAsym_sum_roi(i+input.number_of_baselines).auc_around_mid_point;
    C1(i) = auc_BS_sum_roi(i+input.number_of_baselines).auc_around_mid_point;
end
mean_auc_glucose2_STAsym_sum_roi = mean(A1);
 std_auc_glucose2_STAsym_sum_roi = std(A1);
mean_auc_glucose2_BS_sum_roi     = mean(C1);
 std_auc_glucose2_BS_sum_roi     = std(C1);

%% save all variables in .mat format
clear i ii jj kk
cd(directory)
save( [num2str(input.number_of_scans), 'scans_', ...
       num2str(input.ppm_mid_point_of_auc), 'ppm_', ...
       num2str(input.regularization_weight),'regF_', ...
       num2str(input.Rsqr_threshold),'Rsqr_',datestr(datetime,'ddmmyyyyHHMM'),'.mat'])

%% print everything in txt and other formats
print_script

%% make all plots
plots
% close all
