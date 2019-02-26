%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First approach: z-spectrum of AVERAGE of baselines' raw z-spectra %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (better approach)
%output: mean_baseline_sum_roi_raw                 interpolated_avg_baseline_sum_roi
%        mean_baseline_Rsqr_sum_roi_raw            interpolated_avg_baseline_Rsqr_sum_roi
%        mean_baseline_rois_raw                    interpolated_avg_baseline_rois
%        mean_baseline_Rsqr_rois_raw               interpolated_avg_baseline_Rsqr_rois
%        mean_baseline_pixels_raw                  interpolated_avg_baseline_pixels 
%normalized values used
%%%%%%%%%%%%%%
%for sum_ROI %
%%%%%%%%%%%%%%
% We don't really care what number of cest scan is used in the next line 
% since all scans have the same number of ppm values
baseline = zeros(1,data_in_ppm_sets(1).number_of_ppm_values);
%finds the mean value over the baseline scans at each ppm point
for ii = 1 : input.number_of_baseline_sets %loop over number of baseline measurements
  for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
    baseline(jj) = baseline(jj) + z_spectra_raw(ii).normalized_sum_roi(jj);
  end
end
mean_baseline_sum_roi_raw = baseline / input.number_of_baseline_sets;
clear baseline
% interpolation of average baseline z-spectrum
[avg_baseline_interpolation_sum_roi] = ... %gets normilized in B_spline_and_STAsym function
 B_spline_Rsqr_STAsym(data_in_ppm_sets(1).ppm_values, mean_baseline_sum_roi_raw,...
                    input.regularization_weight, input.auto_or_custom,input.sign_st);
                             
%%%%%%%%%%%%%%%%%%%
%for Rsqr_sum_ROI %
%%%%%%%%%%%%%%%%%%%
% We don't really care what number of cest scan is used in the next line 
% since all scans have the same number of ppm values
baseline = zeros(1,data_in_ppm_sets(1).number_of_ppm_values);
%finds the mean value over the baseline scans at each ppm point
for ii = 1 : input.number_of_baseline_sets %loop over number of baseline measurements
  for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
    baseline(jj) = baseline(jj) + z_spectra_raw(ii).normalized_Rsqr_sum_roi(jj);
  end
end
mean_baseline_Rsqr_sum_roi_raw = baseline / input.number_of_baseline_sets;
clear baseline
% interpolation of average baseline z-spectrum
[avg_baseline_interpolation_Rsqr_sum_roi] = ... %gets normilized in B_spline_and_STAsym function
 B_spline_Rsqr_STAsym(data_in_ppm_sets(1).ppm_values, mean_baseline_Rsqr_sum_roi_raw,...
                    input.regularization_weight, input.auto_or_custom,input.sign_st);                             
                            
%%%%%%%%%%%%%%%%%%%%%%
%for individual ROIs %
%%%%%%%%%%%%%%%%%%%%%%                             
%  z-spectrum of mean of individual ROIs of baselines
if masks.number_of_rois > 1 
 % We don't really care what number of cest scan is used in the next line 
 % since all scans have the same number of ppm values, same in "jj 'for' loop"
  baseline_rois = zeros(masks.number_of_rois,data_in_ppm_sets(1).number_of_ppm_values);
  for kk = 1 : masks.number_of_rois
    for jj = 1 : data_in_ppm_sets(1).number_of_ppm_values
      for ii = 1 : input.number_of_baseline_sets
        baseline_rois(kk,jj) = baseline_rois(kk,jj) + ...
                               z_spectra_raw(ii).normalized_rois(kk,jj);
      end
    end
    mean_baseline_rois_raw = baseline_rois / input.number_of_baseline_sets; 
    %the "mean_baseline_rois_normalized_raw" will be calculated in the
    %"B-spline" function in a while!
    
    % interpolation of average baseline z-spectrum of each individual ROI 
  [avg_baseline_interpolation_rois(kk)] = ... %gets normilized in B_spline_and_STAsym function
   B_spline_Rsqr_STAsym(data_in_ppm_sets(1).ppm_values, mean_baseline_rois_raw(kk,:),...
                      input.regularization_weight, input.auto_or_custom,input.sign_st); 
  end
  clear baseline_rois
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  %for Rsqr individual ROIs %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  baseline_Rsqr_rois = zeros(masks.number_of_rois,cest_scan(1).number_of_ppm_values);
  for kk = 1 : masks.number_of_rois
    for jj = 1 : data_in_ppm_sets(1).number_of_ppm_values
      for ii = 1 : input.number_of_baseline_sets
        baseline_Rsqr_rois(kk,jj) = baseline_Rsqr_rois(kk,jj) + ...
                                  z_spectra_raw(ii).normalized_Rsqr_rois(kk,jj);
      end
    end
    mean_baseline_Rsqr_rois_raw = baseline_Rsqr_rois / input.number_of_baseline_sets; 
    %the "mean_baseline_rois_normalized_raw" will be calculated in the
    %"B-spline" function in a while!
  end
    % interpolation of average baseline z-spectrum of each individual ROI 
  [avg_baseline_interpolation_Rsqr_rois(kk)] = ... %gets normilized in B_spline_and_STAsym function
   B_spline_Rsqr_STAsym(data_in_ppm_sets(1).ppm_values, mean_baseline_Rsqr_rois_raw(kk,:),...
                      input.regularization_weight, input.auto_or_custom,input.sign_st);  
  clear baseline_Rsqr_rois
end

%%%%%%%%%%%%%%%%%%%%%%%
% for pixel wise data %
%%%%%%%%%%%%%%%%%%%%%%%
baseline = zeros(data_in_ppm_sets(1).number_of_ppm_values,size(interpolated_pixels,2));
g = waitbar(0,'Baseline averaging on pixel data...');
%finds the mean value over the baseline scans at each ppm point
for ll = 1 : size(interpolated_pixels,2) 
  waitbar(ll/size(interpolated_pixels,2))  
  for ii = 1 : input.number_of_baseline_sets %loop over number of baseline measurements
    for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
      baseline(jj,ll) = baseline(jj,ll) + interpolated_pixels(ii,ll).non_interp_normalized_data(jj);
    end
  end

mean_baseline_pixels_raw = baseline / input.number_of_baseline_sets;
% interpolation of average baseline z-spectrum
[avg_baseline_interpolation_pixels(ll)] = ... %gets normilized in B_spline_and_STAsym function
 B_spline_Rsqr_STAsym(data_in_ppm_sets(1).ppm_values, mean_baseline_pixels_raw(:,ll),...
                    input.regularization_weight, input.auto_or_custom,input.sign_st);
end
close (g)
clear baseline

disp('avg_of_raw_baselines_and_interpolation DONE');