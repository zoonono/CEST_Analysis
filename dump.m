%% mask image with roi

nnz_masked_sets = nnz(masks.sum_roi_resized);

for ii = 1 : length(data_in_ppm_sets)
  for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
    %temporarely created variable image_temp for the multiplication    
    image_temp(:,:) = data_in_ppm_sets(ii).image_in_gray_scale(:,:,jj);
    %scans masked only with RESIZED sum_ROI mask
    masked_sets(ii).image_in_gray_scale(:,:,jj) = immultiply(image_temp, masks.sum_roi_resized); 
    
    raw_spectra_roi(ii,jj) = sum(sum(masked_sets(ii).image_in_gray_scale(:,:,jj)))/nnz_masked_sets;
  end
end





% % %% put all ppm values in a timeline stack of images
% % i=1;
% % 
% % for ii = 1: length(data_in_ppm_sets)
% %   for jj = 1: data_in_ppm_sets(ii).number_of_ppm_values  
% %     images{i} = data_in_ppm_sets(ii).image_in_gray_scale(:,:,jj);
% %     i = i + 1;
% %   end
% % end
% % 
% % Imatrix = []; 
% % for i = 1:length(data_in_ppm_sets)* data_in_ppm_sets(ii).number_of_ppm_values
% %     Imatrix = cat(3, Imatrix, images{i});
% % end
% % 
% % for jj = 1: data_in_ppm_sets(ii).number_of_ppm_values  
% %     figure
% %     imshow(data_in_ppm_sets(ii).image_in_gray_scale(:,:,jj))
% % end
% % 
