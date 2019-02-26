%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Afroditi Eleftheriou                                                    %
% 22.12.2017                                                              %
%                                                                         %
% Print all data in txt files using the function 'print_in_txt'           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print RAW data (mean values of acquired ppm values) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first load data from struct to temporally used vectors "sum_roi" 
%and "normalized_sum_roi" 
for ii = 1 : length(data_in_ppm_sets)
  for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
    if masks.number_of_rois > 1 
      for kk = 1 : masks.number_of_rois  
        rois(kk,ii,jj) = z_spectra_raw(ii).rois(kk,jj);
        normalized_rois(kk,ii,jj) = ...
                   z_spectra_raw(ii).normalized_rois(kk,jj);
      end
    end
    sum_roi(ii,jj) = z_spectra_raw(ii).mean_sum_roi(jj);
    normalized_sum_roi(ii,jj) = z_spectra_raw(ii).normalized_sum_roi(jj);    
  end
end

% path = strcat (uigetdir, '\');

if masks.number_of_rois > 1 
  for kk = 1 : masks.number_of_rois
    filename = strcat('raw_mean_of_roi_No',num2str(kk),'.txt');
    %"reshape" makes a cell of e.g. 1x20x47 to a vector 20x47
    print_in_txt(data_in_ppm_sets(1).ppm_values, length(data_in_ppm_sets), ...
                 reshape(rois(kk,:,:), ... 
                         [length(data_in_ppm_sets), data_in_ppm_sets(1).number_of_ppm_values]), ...
                 filename, path);
    filename = strcat('raw_normalized_of_roi_No',num2str(kk),'.txt');
    print_in_txt(data_in_ppm_sets(1).ppm_values,length(data_in_ppm_sets), ...
             reshape(normalized_rois(kk,:,:), ...
                     [length(data_in_ppm_sets), data_in_ppm_sets(1).number_of_ppm_values]), ...
             filename, path);
  end
end

% filename = 'raw_sum_roi.txt';
% print_in_txt(cest_scan(1).ppm_values, input.number_of_scans, sum_roi, filename, path);        
         
% filename = 'raw_normalized_sum_roi.txt';
% print_in_txt(cest_scan(1).ppm_values,input.number_of_scans, ...
%              normalized_sum_roi, filename, path);
         
clear rois   normalized_rois   sum_roi normalized_sum_roi
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print INTERPOLATED data (mean values of acquired ppm values) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first load data from struct to temporally used vectors "sum_roi" 
%and "normalized_sum_roi" 
for ii = 1 : length(data_in_ppm_sets)
  for jj = 1 : length(interpolated_sum_roi(ii).ppm)
    if masks.number_of_rois > 1 
      for kk = 1 : masks.number_of_rois  
        normalized_rois(kk,ii,jj) = ...
                   interpolated_rois(kk,ii).z_spectrum(jj);
      end
    end
    normalized_sum_roi(ii,jj) = interpolated_sum_roi(ii).z_spectrum(jj);  
  end
end

if masks.number_of_rois > 1 
  for kk = 1 : masks.number_of_rois
    %"reshape" makes a cell of e.g. 1x20x47 to a vector 20x47
    filename = strcat('interpolated_normalized_of_roi_No',num2str(kk),'.txt');
    print_in_txt(interpolated_sum_roi(1).ppm,length(data_in_ppm_sets), ...
             reshape(normalized_rois(kk,:,:), ...
                     [length(data_in_ppm_sets), length(interpolated_sum_roi(1).ppm)]), ...
             filename, path);
  end
end         
         
% filename = 'interpolated_normalized_sum_roi.txt';
% print_in_txt(interpolated_sum_roi(1).ppm,input.number_of_scans, ...
%              normalized_sum_roi, filename, path);
       
%datestr(datetime,'ddmmyyyyHHMMSS')