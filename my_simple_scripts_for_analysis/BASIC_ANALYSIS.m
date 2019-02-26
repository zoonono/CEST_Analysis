read_files_and_rois

%% mask image with roi

%choose manually with mask to use
mask = masks.sum_roi_resized;
%non zero pixels of this mask, for a correct calculation of the mean
nnz_masked_sets = nnz(mask);

for ii = 1 : length(data_in_ppm_sets)
  for jj = 1 : data_in_ppm_sets(ii).number_of_ppm_values
    %temporarely created variable image_temp for the multiplication    
    image_temp(:,:) = data_in_ppm_sets(ii).image_in_gray_scale(:,:,jj);
    %scans masked only with RESIZED sum_ROI mask
    masked_sets(ii).image_in_gray_scale(:,:,jj) = immultiply(image_temp, mask); 
    
    raw_spectra_roi(ii,jj) = sum(sum(masked_sets(ii).image_in_gray_scale(:,:,jj)))/nnz_masked_sets;
  end
end

%% plot timeline of 1.2ppm normalized to first baseline

f = figure;
set(f, 'PaperOrientation', 'portrait'); %'PaperSize', [21 29.7],

f1 = subplot(2,1,1);
signal1 = (raw_spectra_roi(:,3))/raw_spectra_roi(1,3);
plot(experiment_time(:),signal1)
min1 = min(signal1);
max1 = max(signal1)+ 0.05*max(signal1);
   line([10,11], [min1 + 0.05*min1,min1 + 0.05*min1],'Color', 'black','LineWidth', 1.5) %([x,x],[y,y])
   line([25,70], [min1 + 0.05*min1,min1 + 0.05*min1],'Color', 'black','LineWidth', 1.5) %([x,x],[y,y])
   line([83,130],[min1 + 0.05*min1,min1 + 0.05*min1],'Color', 'black','LineWidth', 1.5) %([x,x],[y,y])
   line([10,10],   [min1,max1],'Color', [0.8 0.8 0.8],'LineWidth', 0.5) %([x,x],[y,y])
   line([110,110], [min1,max1],'Color', [0.8 0.8 0.8],'LineWidth', 0.5) %([x,x],[y,y])
   line([120,120], [min1,max1],'Color', [0.8 0.8 0.8],'LineWidth', 0.5) %([x,x],[y,y])
   line([130,130], [min1,max1],'Color', [0.8 0.8 0.8],'LineWidth', 0.5) %([x,x],[y,y])
   line([140,140], [min1,max1],'Color', [0.8 0.8 0.8],'LineWidth', 0.5) %([x,x],[y,y])
   ylim ([min1 max1]) 
%plot(experiment_time(21:220),raw_spectra_roi(21:220,3)/raw_spectra_roi(21,3))
ax = ancestor(f1, 'axes');
   ax.YAxis.FontSize = 14;
   ax.XAxis.FontSize = 14;

   x3 = xlabel('Time (min)');
   x3.FontSize = 13;

   y3 = ylabel('glucoCEST intensity');
   y3.FontSize = 13;

   t = title('glucoCEST intensity normalized to baseline, @1.2ppm');
   t.FontSize = 14;
   
   
   
f2 = subplot(2,1,2);
signal2 = 100*(raw_spectra_roi(:,3)-raw_spectra_roi(1,3))/raw_spectra_roi(1,3);
plot(experiment_time(:),signal2)
min2 = min(signal2);
max2 = max(signal2)+ 0.05*max(signal2);
   line([10,11], [min2 + 0.05*min2,min2 + 0.05*min2],'Color', 'black','LineWidth', 1.5) %([x,x],[y,y])
   line([25,70], [min2 + 0.05*min2,min2 + 0.05*min2],'Color', 'black','LineWidth', 1.5) %([x,x],[y,y])
   line([83,130],[min2 + 0.05*min2,min2 + 0.05*min2],'Color', 'black','LineWidth', 1.5) %([x,x],[y,y])
   line([10,10],   [min2,max2],'Color', [0.8 0.8 0.8],'LineWidth', 0.5) %([x,x],[y,y])
   line([110,110], [min2,max2],'Color', [0.8 0.8 0.8],'LineWidth', 0.5) %([x,x],[y,y])
   line([120,120], [min2,max2],'Color', [0.8 0.8 0.8],'LineWidth', 0.5) %([x,x],[y,y])
   line([130,130], [min2,max2],'Color', [0.8 0.8 0.8],'LineWidth', 0.5) %([x,x],[y,y])
   line([140,140], [min2,max2],'Color', [0.8 0.8 0.8],'LineWidth', 0.5) %([x,x],[y,y])
   ylim ([min2 max2])
%plot(experiment_time(21:220),raw_spectra_roi(21:220,3)/raw_spectra_roi(21,3))
ax2 = ancestor(f2, 'axes');
   ax2.YAxis.FontSize = 14;
   ax2.XAxis.FontSize = 14;

   x2 = xlabel('Time (min)');
   x2.FontSize = 13;

   y2 = ylabel('Baseline Subtraction (%)');
   y2.FontSize = 13;

   t = title('glucoCEST Baseline Subtraction, normalized to baseline @1.2ppm');
   t.FontSize = 14;   
   
   cd('D:\DATA_glucoCEST\omg_ins_glc_gal');
   print('-bestfit', [input.mouse_name,'.pdf'],'-dpdf')
