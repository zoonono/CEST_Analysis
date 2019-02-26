%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%-%         P I X E L - W I S E   I N T E R P O L A T I O N  ! !        %-%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%interpolate pixel-wise z-spectra (ONLY FOR PIXEL WITHIN sum_roi !!!) 
%Also calculates R^2 of each pixel ("goodness of interpolation")
h = waitbar(0,'Pixel-wise z-spectra interpolation...');
for ii = 1 : length(data_in_ppm_sets)
  waitbar(ii/length(data_in_ppm_sets))
  number_of_pixel = 0;  %within each scan it'll calculate the number of non zero pixels
  for jj = 1 : size(masked_scan(ii).sum_roi,1)
    for kk = 1 : size(masked_scan(ii).sum_roi,2)
      if masks.sum_roi_resized(jj,kk) > 0 %check if the pixel is within the sum_ROI
        number_of_pixel = number_of_pixel + 1;
        %temporary struct name "interpolated", will be named "interpolated_pixels"  
        %the data get NORMALIZED in the "B_spline_Rsqr_STAsym" function
        [interpolated] = B_spline_Rsqr_STAsym(data_in_ppm_sets(ii).ppm_values, ... 
                         reshape(masked_scan(ii).sum_roi(jj,kk,:), ...
                         [1,data_in_ppm_sets(ii).number_of_ppm_values]),...
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