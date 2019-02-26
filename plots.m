%make some values in the correct format
date = [cest_scan(1).data_date_and_time(9:10),'.',...
        cest_scan(1).data_date_and_time(6:7),'.',...
        cest_scan(1).data_date_and_time(1:4)];
minutes = num2str(scan_duration_min_sec(1)); % temporal resolution "min" part
seconds = num2str(scan_duration_min_sec(2)); % temporal resolution "sec" part

%% create figure
figure_1 = figure; %opening an empty figure space where all graphs and textboxes will go
set(figure_1, 'PaperSize', [21 29.7], 'PaperOrientation', 'landscape');
%first textbox (Date etc)
h1 = subplot(3,4,1); %3lines, 3 columns, position 1
h1.Position = [0.02 0.82 0.2 0.1]; %[left bottom width height]
axis off 
t1= text;
t1.FontSize = 11 ; 
t1.Interpreter = 'tex'; 
t1.String = { date,...
             [' mouse: ',input.mouse_name],...
             input.administration,...
             ['fasted: ', input.fasted],...
             ['# baseline scans: ', num2str(input.number_of_baseline_sets)],...
             ['# post admin. scans: ', num2str(input.number_of_post_inj_sets)],...
             ['# ppm values: ', num2str(data_in_ppm_sets(1).number_of_ppm_values)],...
             ['temp. resolution: ', minutes, 'min ', seconds,'sec'],...
             ' ',...
             ['ST-Asym: mean \pm std = ', ... 
               num2str(mean_auc_glucose_STAsym_sum_roi),' \pm ', ...
               num2str(std_auc_glucose_STAsym_sum_roi),...
              ' (',num2str(100*std_auc_glucose_STAsym_sum_roi/mean_auc_glucose_STAsym_sum_roi),'%)'],...
             ['BS-method: mean \pm std = ', ...
               num2str(mean_auc_glucose_BS_sum_roi),' \pm ', ...
               num2str(std_auc_glucose_BS_sum_roi),...
              ' (',num2str(100*std_auc_glucose_BS_sum_roi/mean_auc_glucose_BS_sum_roi),'%)']};
 
%%
%second textbox (pre-sat. parameters etc)
h2 = subplot(3,4,2); %3lines, 3 columns, position 2
h2.Position = [0.23 0.857 0.2 0.1]; %[left bottom width height]
axis off
t2=text;      
t2.FontSize = 11 ;
t2.Interpreter = 'tex';                     
t2.String = {['pre-sat. strength: ',num2str(cest_scan(1).B1_power),'\muT'],...
             ['pre-sat. duration: ',num2str(cest_scan(1).saturation_pulse_duration/1000),'s'],...         
             ['TE = ',num2str(cest_scan(1).TE),'ms'],...
             ['TR = ',num2str(cest_scan(1).TR),'ms'],...
             ['RARE = ',num2str(cest_scan(1).RARE)],...
             ['slice thickness = ',num2str(cest_scan(1).slice_thickness),'mm'],...
             ['FOV = ', num2str(cest_scan(1).FOV(1,1)),'mm x ',num2str(cest_scan(1).FOV(1,2)),'mm'],...
             ['matrix = ', num2str(cest_scan(1).size_x),' x ', num2str(cest_scan(1).size_y)]};        

%% plot ppm VS interpolated normalized mean of sum ROI
h34 = subplot(3,4,[3 4]);
h34.Position = [0.55 0.75 0.45 0.22]; %[left bottom width height]
%f2 = figure('Name','Interpolated z-spectra of sum ROI','NumberTitle','off');
for i = 1:length(data_in_ppm_sets)
    plot(interpolated_sum_roi(i).ppm, interpolated_sum_roi(i).z_spectrum,'-');
    legendInfo{i}=['scan ',num2str(i)];
    hold on 
end
plot(interpolated_sum_roi(1).ppm, avg_baseline_interpolation_sum_roi.z_spectrum, '--r')
title('Interpolated z-spectra of ROI')
xlabel('Saturation frequency (p.p.m.)')
ylabel('Nomalized Intensity')
xlim([min(interpolated_sum_roi(1).ppm) max(interpolated_sum_roi(1).ppm)])
ylim([min(interpolated_sum_roi(1).z_spectrum) max(interpolated_sum_roi(1).z_spectrum)])

% legendInfo{i+1} = 'avg baseline';
% lgd = legend('show');
% lgd.String = legendInfo;
% lgd.Location = 'southwest';
% lgd.Box = 'off';

clear legendInfo

%% plot interpolated z-spectrum ST_ASYMMETRY CURVE of mean of sum ROI
h56 = subplot(3,4,[5 6]);
h56.Position = [0.02 0.40 0.45 0.22]; %[left bottom width height]
%f3 = figure('Name','Interpolated ST asymmetry curve of sum ROI');
for i = 1 : length(data_in_ppm_sets)
    plot(interpolated_sum_roi(i).ppm_positive, interpolated_sum_roi(i).STAsym,'-');
    legendInfo{i}=['scan ',num2str(i)];
    hold on
end
title('Interpolated ST Asymmetry curve of ROI')
xlabel('Saturation frequency (p.p.m.)')
ylabel('Intensity (arbitrary units)')
% lgd = legend('show');
% lgd.String = legendInfo;
% lgd.Box = 'off';

clear legendInfo
hold off

%% plot AUC of ST Asymmetry curve analysis
h78 = subplot(3,4,[7 8]);
h78.Position = [0.55 0.40 0.45 0.22]; %[left bottom width height]
%f8 = figure('Name','Area under the curve of ST Asymmetry curve');
for i_time = 1: length(data_in_ppm_sets)
    y(i_time) = auc_STAsym_sum_roi(i_time).auc_around_mid_point/ ...
                mean_auc_glucose_STAsym_sum_roi - 1;
end
plot(experiment_time, y,'-o');
title('AUC ST Asymmetry curve')
xlabel('time (sec)')
ylabel('Normalized Intensity')

%% plot z-spectra of interpolated-mean-baseline SUBTRACTION of sum ROI
h910 = subplot(3,4, [9 10]);
h910.Position = [0.02 0.05 0.45 0.22]; %[left bottom width height]
%f5 = figure('Name','z-spectra of interpolated-mean-baseline SUBTRACTION of sum ROI');
for i = 1 : length(data_in_ppm_sets)
   plot(interpolated_sum_roi(1).ppm, BS_interpolated_mean_sum_roi(i,:),'-');
   legendInfo{i} = ['scan ',num2str(i)];
   hold on
end
title('z-spectra of interpolated-mean-baseline SUBTRACTION of sum ROI')
xlabel('Saturation frequency (p.p.m.)')
ylabel('Intensity (arbitrary units)')

% lgd = legend('show');
% lgd.String = legendInfo;
% lgd.Box = 'off';

clear legendInfo
hold off
%% plot AUC of baseline subtraction analysis
h1112 = subplot(3,4,[11 12]);
h1112.Position = [0.55 0.05 0.45 0.22]; %[left bottom width height]
%f7 = figure('Name','Area under the curve of average baseline subtraction');
for i_time = 1: length(data_in_ppm_sets)
    y(i_time) = auc_BS_sum_roi(i_time).auc_around_mid_point/ ...
                mean_auc_glucose_BS_sum_roi -1;
end
% int  = input.number_of_baselines + 1;
% last = int + input.number_of_post_inj_scans -1;
%plot(experiment_time(int:last), y(int:last),'-o');
plot(experiment_time, y,'-o');
title('AUC average baseline subtraction')
xlabel('time (sec)')
ylabel('Normalized Intensity')

clear y

%% save figure with all needed graphs
print('-bestfit',[num2str(length(data_in_ppm_sets)), 'scans_', ...
       num2str(input.ppm_mid_point_of_auc), 'ppm_', ...
       num2str(input.regularization_weight),'regF_', ...
       num2str(input.Rsqr_threshold),'Rsqr_',datestr(datetime,'ddmmyyyyHHMM'),'.pdf'],'-dpdf')
        

% % %% images
% % % R squared image
% % for i = 1: input.number_of_scans
% % Rsqr_figure = squeeze(image_Rsqr(i,:,:));
% % figure('Name',[ 'R squared map of scan ', num2str(i)],'NumberTitle','off');
% % imshow(Rsqr_figure, [min(Rsqr_figure(Rsqr_figure>0)) 1]);
% % end
% % 
% % % area under the curve images (glucose contrast images)
% % % for "auc_STAsym_pixels"
% % 
% % % image_AUC_around_glc_peak = zeros(ii,cest_scan(ii).size_x,cest_scan(ii).size_y);
% % % image_AUC_zero_to_peak = zeros(ii,cest_scan(ii).size_x,cest_scan(ii).size_y);
% % for ii = 1 : input.number_of_scans
% %   for ll = 1 : size(interpolated_pixels,2)
% %     x = interpolated_pixels(ii,ll).pixel_position(1);
% %     y = interpolated_pixels(ii,ll).pixel_position(2);
% %     image_AUC_around_glc_peak(ii,x,y) = auc_STAsym_pixels(ii,ll).auc_around_mid_point;
% %     image_AUC_zero_to_peak(ii,x,y) = auc_STAsym_pixels(ii,ll).auc_from_zero;
% %   end
% % end
% % 
% % %image AUC around glucose peak
% % for ii = 1 : input.number_of_scans
% %     glc_image = squeeze(image_AUC_around_glc_peak(ii,:,:));
% %     glc_image_resized = imresize(glc_image, ...
% %                 [size(anatomical_scan.image_in_gray_scale,1) size(anatomical_scan.image_in_gray_scale,2)]);
% %     figure('Name',['glucose contrast image, AUC around glucose peak, scan ', num2str(ii)],...
% %            'NumberTitle','off');
% %     imshow(anatomical_scan.image_in_gray_scale);
% %     colormap gray
% %     freezeColors
% %     hold on
% %     h = imshow((medfilt2(double(glc_image_resized),[2 2])), [0 0.15]);
% %     colormap(gca, 'hot')
% %     alpha 0.7
% %     colorbar
% %     hold off
% % end    
% % 
% % %image AUC from zero
% % for ii = 1 : input.number_of_scans
% %     glc_image = squeeze(image_AUC_zero_to_peak(ii,:,:));
% %     glc_image_resized = imresize(glc_image, ...
% %                 [size(anatomical_scan.image_in_gray_scale,1) size(anatomical_scan.image_in_gray_scale,2)]);
% %     figure('Name',['glucose contrast image, AUC from zero, scan ', num2str(ii)],...
% %            'NumberTitle','off');
% %     imshow(anatomical_scan.image_in_gray_scale);
% %     colormap gray
% %     freezeColors
% %     hold on
% %     %B = medfilt2(A) performs median filtering of the matrix A in two dimensions
% %     imshow((medfilt2(double(glc_image_resized),[2 2])) , [0 0.3]) 
% %     colormap(gca, 'hot')
% %     alpha 0.7
% %     colorbar
% %     hold off
% % end    

% %%graphs that are not thaaaaat needed!

% %% plot z-spectra of mean-of-interpolated-baselines SUBTRACTION of sum ROI
% f6 = figure('Name','z-spectra of mean-of-interpolated-baselines SUBTRACTION of sum ROI');
% for i = 1 : input.number_of_scans
%    plot(interpolated_sum_roi(1).ppm, BS_mean_of_interpolated_sum_roi(i,:),'-');
%    legendInfo{i} = ['scan ',num2str(i)];
%    hold on
% end
% title('z-spectra of mean-of-interpolated-baselines SUBTRACTION of ROI')
% xlabel('Saturation frequency (p.p.m.)')
% ylabel('Intensity (arbitrary units)')
% 
% lgd = legend('show');
% lgd.String = legendInfo;
% lgd.Box = 'off';
% 
% clear legendInfo
% hold off
% 
% %% plot ppm VS normalized mean of sum ROI
% f1 = figure('Name','Raw z-spectra of sum ROI','NumberTitle','off'); 
% for i = 1:input.number_of_scans
%     plot(cest_scan(i).ppm_values, z_spectra_raw(i).normalized_sum_roi,'-');
%     legendInfo{i}=['scan ',num2str(i)];
%     hold on 
% end
% plot(cest_scan(1).ppm_values,mean_baseline_sum_roi_raw, '-*r')
% title('Raw z-spectra of ROI')
% xlabel('Saturation frequency (p.p.m.)')
% ylabel('Intensity (arbitrary units)')
% xlim([min(cest_scan(1).ppm_values) max(cest_scan(1).ppm_values)])
% 
% legendInfo{i+1} = 'avg baseline';
% lgd = legend('show');
% lgd.String = legendInfo;
% lgd.Location = 'southwest';
% lgd.Box = 'off';
% 
% 
% clear legendInfo
% hold off
% %% plot RAW z-spectrum SUBTRACTION average baseline, for sum ROI
% f4 = figure('Name','Average baseline subtraction, raw z-spectra of sum ROI','NumberTitle','off');
% for i = 1 : input.number_of_scans
%    plot(cest_scan(i).ppm_values, BS_RAW_sum_roi(i,:),'-');
%    legendInfo{i} = ['scan ',num2str(i)];
%    hold on
% end
% 
% title('Raw z-spectra subtraction (avg baseline) of ROI')
% xlabel('Saturation frequency (p.p.m.)')
% ylabel('BS intensity (arbitrary units)')
% 
% lgd = legend('show');
% lgd.String = legendInfo;
% lgd.Box = 'off';
% 
% clear legendInfo
% hold off

