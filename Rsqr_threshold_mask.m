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

for ii = 1 : length(data_in_ppm_sets)
  for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
    %temporarely created variable image_2dseq for the multiplication    
    image_2dseq (:,:) = data_in_ppm_sets(ii).image_2dseq(:,:,jj);

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