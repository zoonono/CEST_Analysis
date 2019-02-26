%normalized values used
%output: BS_RAW_sum_roi
%        BS_interpolated_mean_Rsqr_sum_roi
%        BS_interpolated_mean_sum_roi
%        BS_RAW_rois
%        BS_interpolated_mean_Rsqr_rois
%        BS_RAW_pixels
%        BS_interpolated_mean_pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for sum_ROI and Rsqr_sum_ROI %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%RAW baseline subtraction (BS) (in RAW data)
BS_RAW_sum_roi = zeros(length(data_in_ppm_sets),data_in_ppm_sets(1).number_of_ppm_values); 
for ii = 1 : length(data_in_ppm_sets)
 for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
  %mean_baseline_sum_roi_raw has normalized values
  BS_RAW_sum_roi(ii,jj) = z_spectra_raw(ii).normalized_sum_roi(jj) - ...
                          mean_baseline_sum_roi_raw(jj);
 end
end

%Rsqr_sum_roi
%interpolated mean-of-raw-baselines subtraction (1rst approach)
BS_interpolated_mean_Rsqr_sum_roi = zeros(length(data_in_ppm_sets),length(interpolated_sum_roi(1).ppm));
for ii = 1 : length(data_in_ppm_sets)
 for jj = 1 : length(interpolated_sum_roi(1).ppm)
  BS_interpolated_mean_Rsqr_sum_roi(ii,jj) = interpolated_Rsqr_sum_roi(ii).z_spectrum(jj) - ...
                                             avg_baseline_interpolation_Rsqr_sum_roi.z_spectrum(jj);
 end
end
%sum_roi
%interpolated mean-of-raw-baselines subtraction (1rst approach)
BS_interpolated_mean_sum_roi = zeros(length(data_in_ppm_sets),length(interpolated_sum_roi(1).ppm));
for ii = 1 : length(data_in_ppm_sets)
 for jj = 1 : length(interpolated_sum_roi(1).ppm)
  BS_interpolated_mean_sum_roi(ii,jj) = interpolated_sum_roi(ii).z_spectrum(jj) - ...
                                        avg_baseline_interpolation_sum_roi.z_spectrum(jj);
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for individual ROIs and Rsqr_individual_ROIs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if masks.number_of_rois > 1 
BS_RAW_rois = zeros(masks.number_of_rois, length(data_in_ppm_sets), data_in_ppm_sets(1).number_of_ppm_values);
BS_interpolated_mean_Rsqr_rois = zeros(masks.number_of_rois,length(data_in_ppm_sets),length(interpolated_sum_roi(1).ppm));
% BS_mean_of_interpolated_Rsqr_rois = zeros(masks.number_of_rois,input.number_of_scans,length(interpolated_sum_roi(1).ppm));
for kk = 1 : masks.number_of_rois
 
  %RAW baseline subtraction (BS)  
  for ii = 1 : length(data_in_ppm_sets)
   for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
    BS_RAW_rois(kk,ii,jj) = z_spectra_raw(ii).normalized_rois(kk,jj) - mean_baseline_rois_raw(kk,jj);
   end
  end
  
  %Rsqr_rois
  %interpolated mean-of-raw-baselines subtraction (1rst approach)
  for ii = 1 : length(data_in_ppm_sets)
   for jj = 1 : length(interpolated_sum_roi(1).ppm)
     BS_interpolated_mean_Rsqr_rois(kk,ii,jj) = interpolated_Rsqr_rois(kk,ii).z_spectrum(jj) - ...
                                                avg_baseline_interpolation_Rsqr_rois(kk).z_spectrum(jj);
   end
  end
    %rois
  %interpolated mean-of-raw-baselines subtraction (1rst approach)
  for ii = 1 : length(data_in_ppm_sets)
   for jj = 1 : length(interpolated_sum_roi(1).ppm)
     BS_interpolated_mean_rois(kk,ii,jj) = interpolated_rois(kk,ii).z_spectrum(jj) - ...
                                           avg_baseline_interpolation_rois(kk).z_spectrum(jj);
   end
  end
end
end

%%%%%%%%%%%%%%
% for pixels %
%%%%%%%%%%%%%%
BS_RAW_pixels = zeros(length(data_in_ppm_sets), ...
                      data_in_ppm_sets(1).number_of_ppm_values,...
                      size(interpolated_pixels,2));
BS_interpolated_mean_pixels = zeros(length(data_in_ppm_sets),...
                                    length(interpolated_sum_roi(1).ppm),...
                                    size(interpolated_pixels,2));
h = waitbar(0,'Baseline subtraction on pixel data...');                                
for ll = 1 : size(interpolated_pixels,2) 
  waitbar(ll/size(interpolated_pixels,2))  
  
  %RAW baseline subtraction (BS)  
  for ii = 1 : length(data_in_ppm_sets)
   for jj = 1 : data_in_ppm_sets(1).number_of_ppm_values
    BS_RAW_pixels(ii,jj,ll) = interpolated_pixels(ii,ll).non_interp_normalized_data(jj) - ...
                              mean_baseline_pixels_raw(jj,ll);
   end
  end
  
  %interpolated mean-of-raw-baselines subtraction (1rst approach)
  for ii = 1 : length(data_in_ppm_sets)
   for jj = 1 : length(interpolated_sum_roi(1).ppm)
     BS_interpolated_mean_pixels(ii,jj,ll) = interpolated_pixels(ii,ll).z_spectrum(jj) - ...
                                             avg_baseline_interpolation_pixels(ll).z_spectrum(jj);
   end
  end
end
close(h)