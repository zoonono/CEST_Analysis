%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%-%                                                                     %-%
%-%   M E A N   V A L U E S   A  N D   I N T E R P O L A T I O N  ! !   %-%
%-%                                                                     %-%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%%%%%%%%%%%%%%%
% for sum_ROI %
%%%%%%%%%%%%%%%
%mean value per ppm image, over sum_ROI masked pixel values
for ii = 1 : length(data_in_ppm_sets)
    for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
        z_spectra_raw(ii).mean_sum_roi(jj)= mean2(masked_scan(ii).sum_roi(:,:,jj));
    end
end
%calculate the normalized_mean_sum_roi
for ii = 1 : length(data_in_ppm_sets)
  max_value_of_roi = max (z_spectra_raw(ii).mean_sum_roi(jj));
  for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
  z_spectra_raw(ii).normalized_sum_roi(jj) = ...
                  z_spectra_raw(ii).mean_sum_roi(jj)/ max_value_of_roi;           
  end
end
%INTERPOLATE z_spectra_raw.mean_sum_roi
for ii = 1: length(data_in_ppm_sets)
  [interpolated_sum_roi(ii)] = ...  %gets normilized in B_spline_Rsqr_STAsym function
  B_spline_Rsqr_STAsym(data_in_ppm_sets(ii).ppm_values, z_spectra_raw(ii).mean_sum_roi,...
           input.regularization_weight, input.auto_or_custom,input.sign_st); 
end

%%%%%%%%%%%%%%%%%%%%
% for Rsqr_sum_ROI %
%%%%%%%%%%%%%%%%%%%%
%mean value per ppm image, over Rsqr_sum_ROI masked pixel values
for ii = 1 : length(data_in_ppm_sets)
    for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
        z_spectra_raw(ii).mean_Rsqr_sum_roi(jj)= mean2(masked_scan(ii).Rsqr_sum_roi(:,:,jj));
    end
end
%calculate the normalized_mean_Rsqr_sum_roi
for ii = 1 : length(data_in_ppm_sets)
  max_value_of_roi = max (z_spectra_raw(ii).mean_Rsqr_sum_roi(jj));
  for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
  z_spectra_raw(ii).normalized_Rsqr_sum_roi(jj) = ...
                    z_spectra_raw(ii).mean_Rsqr_sum_roi(jj)/ max_value_of_roi;           
  end
end
%INTERPOLATE z_spectra_raw.mean_Rsqr_sum_roi
for ii = 1: length(data_in_ppm_sets)
  [interpolated_Rsqr_sum_roi(ii)] = ...  %gets normilized in B_spline_Rsqr_STAsym function
  B_spline_Rsqr_STAsym(data_in_ppm_sets(ii).ppm_values, z_spectra_raw(ii).mean_Rsqr_sum_roi,...
           input.regularization_weight, input.auto_or_custom,input.sign_st); 
end
clear max_value_of_roi

%%%%%%%%%%%%%%%%%%%%%%%
% for individual ROIs %
%%%%%%%%%%%%%%%%%%%%%%%
%mean value per ppm image, over individual ROI masked pixel values
if masks.number_of_rois > 1 
  for ii = 1 : length(data_in_ppm_sets)
    for kk = 1 : masks.number_of_rois
      for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
        z_spectra_raw(ii).rois(kk,jj)= mean2(masked_scan(ii).rois(:,:,jj,kk)); 
      end
    end
  end
  %calculate the normalized_rois
  for ii = 1 : length(data_in_ppm_sets)
    for kk = 1 : masks.number_of_rois
      max_values_of_rois(kk) = max (z_spectra_raw(ii).rois(kk,:));
      for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
      z_spectra_raw(ii).normalized_rois(kk,jj) = ...
                        z_spectra_raw(ii).rois(kk,jj)/ max_values_of_rois(kk);
      end
    end
  end
  %INTERPOLATE z_spectra_raw.rois
  for ii = 1: length(data_in_ppm_sets)
    for kk = 1: masks.number_of_rois
    [interpolated_rois(kk,ii)] = ... %gets normilized in B_spline_Rsqr_STAsym function
    B_spline_Rsqr_STAsym(data_in_ppm_sets(ii).ppm_values, z_spectra_raw(ii).rois(kk,:),...
             input.regularization_weight, input.auto_or_custom,input.sign_st);
    end
  end 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % for individual Rsqr_ROIs %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%% (already within "if #rois > 1" statement)
  %mean value per ppm image, over individual ROI masked pixel values
  for ii = 1 : length(data_in_ppm_sets)
    for kk = 1 : masks.number_of_rois
      for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
        z_spectra_raw(ii).Rsqr_rois(kk,jj)= mean2(masked_scan(ii).Rsqr_rois(:,:,jj,kk)); 
      end
    end
  end
  %calculate the normalized_rois
  for ii = 1 : length(data_in_ppm_sets)
    for kk = 1 : masks.number_of_rois
      max_values_of_rois(kk) = max (z_spectra_raw(ii).Rsqr_rois(kk,:));
      for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
      z_spectra_raw(ii).normalized_Rsqr_rois(kk,jj) = ...
                        z_spectra_raw(ii).Rsqr_rois(kk,jj)/ max_values_of_rois(kk);
      end
    end
  end
  %INTERPOLATE z_spectra_raw.rois
  for ii = 1: length(data_in_ppm_sets)
    for kk = 1: masks.number_of_rois
    [interpolated_Rsqr_rois(kk,ii)] = ... %gets normilized in B_spline_Rsqr_STAsym function
    B_spline_Rsqr_STAsym(data_in_ppm_sets(ii).ppm_values, z_spectra_raw(ii).Rsqr_rois(kk,:),...
             input.regularization_weight, input.auto_or_custom,input.sign_st);
    end
  end
end %end of "if masks.number_of_rois > 1" 
