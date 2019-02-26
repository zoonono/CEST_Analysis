%raw data: data_in_ppm_sets 
%roi mask needed: mask.sum_roi
%masked images: masked_scan.sum_roi, .sum_roi_in_gray_scale, .Rsqr_sum_roi

for ii = 1 : length(data_in_ppm_sets)
    for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
        raw_roi(ii,jj)= ...
            sum(sum((masked_scan(ii).sum_roi_in_gray_scale(:,:,jj)))) / ...
            nnz(masked_scan(ii).sum_roi_in_gray_scale(:,:,3));
    end
end


% for i = 1:length(data_in_ppm_sets)
%     plot(masked_scan(ii).sum_roi_in_gray_scale(:,:,3), experiment_time,'-'); 
% end
% plot(interpolated_sum_roi(1).ppm, avg_baseline_interpolation_sum_roi.z_spectrum, '--r')
% title('Interpolated z-spectra of ROI')
% xlabel('Saturation frequency (p.p.m.)')
% ylabel('Nomalized Intensity')
% xlim([min(interpolated_sum_roi(1).ppm) max(interpolated_sum_roi(1).ppm)])
% ylim([min(interpolated_sum_roi(1).z_spectrum) max(interpolated_sum_roi(1).z_spectrum)])
