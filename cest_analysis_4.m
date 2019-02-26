%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Afroditi Eleftheriou                                          %
% Date:     02.01.2019                                                    %
% Version:  4.1                                                           %
%                                                                         %
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       %
%                 %  Main CEST analysis program!  %                       %
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       %
%                                                                         %
% Analysis of data from ParaVision 6.0.1                                  %
% NOTES: run first "read_files_and_rois.m" or import files                %
%        BS = baseline subtraction (seen in nomenclature)                 %
%        ST = saturation transfer (seen in nomenclature)                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Interpolation of pixel-wise data and R^2 calculation
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%-%         P I X E L - W I S E   I N T E R P O L A T I O N  ! !        %-%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
pixel_wise_interpolation

%% Create Rsqr_threshold MASKs and mask data
Rsqr_threshold_mask

%% mean values and interpolation of MASKED data
mean_values_and_interpolation

%%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%-%                                                                     %-%
%-%          B A S E L I N E   S U B T R A C T I O N   ! ! !            %-%
%-%                                                                     %-%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%  

%% Average baseline z-spectrum
avg_of_raw_baselines_and_interpolation

%% Subtraction of average baseline
baseline_subtraction %(BS)

%% AUC
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%-%                                                                     %-%
%-%          A R E A   U N D E R   T H E   C U R V E   ! ! !            %-%
%-%                                                                     %-%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%done on STAsym interpolated data (interpolated_*.STAsym) and 
%baseline-subtraction interpolated data (BS_interpolated_mean_*)
area_under_the_curve

%% Calculating the mean and std of the baselines and post injection sets
mean_and_std_of_baselines

%% save all variables in .mat format
clear i ii jj kk
cd(directory)
save( [num2str(length(data_in_ppm_sets)), 'scans_', ...
       num2str(input.ppm_mid_point_of_auc), 'ppm_', ...
       num2str(input.regularization_weight),'regF_', ...
       num2str(input.Rsqr_threshold),'Rsqr_',datestr(datetime,'ddmmyyyyHHMM'),'.mat'])

%% print everything in txt and other formats
print_script

%% make all plots
plots
% close all
