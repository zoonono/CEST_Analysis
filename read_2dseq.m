%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Afroditi Eleftheriou                                                    %
% 23.11.2017                                                              %
%                                                                         %
% Reads '2dseq' file                                                      %
%                                                                         %
% Input: folder_number                                                    %
%        directory                                                        %
%        parameters_reco.size_x                                           %
%        parameters_reco.size_y                                           %
%        parameters_reco.reco_bit                                         %
%        parameters_reco.number_of_ppm_values                             %
%        parameters_reco.slope                                            %
%                                                                         %
% Output: parameters.image_2dseq                                          %
%         parameters.image_in_gray_scale                                  %
%                                                                         %
%NOTE: next version: insert question for number of slices and read them!  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% These two commented lines are only for testing
%
% folder_number = 17; 
% directory = 'D:\CEST_data\MRI data\mouse_3OMGCEST\20171031_A_3OMGCEST';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[parameters] = read_2dseq(folder_number,directory,par_reco,par_acqp);

% goes to the directory of the scan folder (the name of the folder is a number)
folder_number = num2str(folder_number); %folder_number is the index that is being looped for all scans in the main program
directory_scan = strcat(directory,'\',folder_number); %directory of the scan
directory_2dseq = strcat(directory_scan,'\pdata\1\');
cd(directory_2dseq)

fid=fopen('2dseq');

%store the actual images in 'parameters.image_2dseq'
for k=1:par_acqp.number_of_ppm_values %Loop over each ppm image of a scan
   for i=1:par_reco.size_x %Loop over x dimension of image
       for j=1:par_reco.size_y %Loop over y dimension of image
           parameters.image_2dseq(i,j,k)=fread(fid,1,par_reco.reco_bit);
       end
   end
end
fclose(fid);

%scales all image values for the correction factor
parameters.image_2dseq=parameters.image_2dseq/par_reco.slope;

%gray scale 2dseq image
parameters.image_in_gray_scale = mat2gray(parameters.image_2dseq(:,:,:));

clear fid directory_2dseq