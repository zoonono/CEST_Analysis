%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Afroditi Eleftheriou                                                    %
% 23.11.2017                                                              %
%                                                                         %
% Reads 'acqp' file                                                       %
%                                                                         %
% Input: folder_number                                                    %
%        directory                                                        %
%                                                                         %
% Output: parameters.scan_directory                                       %
%         parameters.institution                                          %
%         parameters.scanner                                              %     
%         parameters.B0_in_MHz                                            %
%         parameters.number_of_ppm_values                                 %
%         parameters.number_of_slices                                     %
%         parameters.slice_thickness                                      %
%         parameters.TE                                                   %
%         parameters.TR                                                   %
%         parameters.RARE                                                 %
%         parameters.ppm_values                                           %
%         parameters.data_date_and_time                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These two commented lines are only for testing
%
% folder_number = 17; 
% directory = 'D:\CEST_data\MRI data\mouse_3OMGCEST\20171031_A_3OMGCEST';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[parameters] = read_acqp(folder_number,directory)

%sets the directory of the scan folder (the name of the folder is a number)
folder_number = num2str(folder_number); %the index that is being looped for all scans in the main program
directory_scan = strcat(directory,'\',folder_number); %directory of the scan
cd(directory_scan);

parameters.scan_directory = directory_scan;
parameters.folder_number = folder_number;
parameters.number_of_ppm_values = 0;
parameters.ppm_values = [];

file_acqp = textread('acqp','%s','delimiter','=','whitespace',''); %reads the 'acqp' file
    
size_acqp = size(file_acqp); %number of lines of the 'acqp' file

%runs through the lines of the 'acqp' file to extract and store some of the info
for i=1:size_acqp(1) 
   if strcmp(file_acqp(i),'##$ACQ_institution')
      parameters.institution=(char(file_acqp(i+2)));   %institution
   end
   if strcmp(file_acqp(i),'##$ACQ_station')            %scanner
      parameters.scanner=(char(file_acqp(i+2)));      
   end
   if strcmp(file_acqp(i),'##$BF1')                    
      % The gyromagnetic ratio: gamma/(2*pi) is:
      % 42.567 MHz/T for Hydrogen-1 
      % 10.705 MHz for Carbon-13
      % 17.235 MHz/T for Phosphorus -31
      % The frequency of a 7T B0 magnetic field is 300.325 MHz. 
      % The reference frequency of CEST is the water frequency, so
      % B0_in_MHz.
      parameters.B0_in_MHz=str2double(char(file_acqp(i+1)));      %B0 in MHz
   end
   if strcmp(file_acqp(i),'##$ACQ_O2_list_size')
      parameters.number_of_ppm_values=round(str2double(char(file_acqp(i+1)))); %number of ppm values
   end
   if strcmp(file_acqp(i),'##$NSLICES')
       parameters.number_of_slices=str2num(char(file_acqp(i+1))); %number of slices
   end
   if strcmp(file_acqp(i),'##$ACQ_slice_thick')
       parameters.slice_thickness=str2num(char(file_acqp(i+1))); %slice thickness
   end   
   if strcmp(file_acqp(i),'##$ACQ_echo_time')
       parameters.TE=str2num(char(file_acqp(i+2))); %echo time in ms
   end 
   if strcmp(file_acqp(i),'##$ACQ_repetition_time')
       parameters.TR=str2num(char(file_acqp(i+2))); %repetition time in ms
   end    
   if strcmp(file_acqp(i),'##$ACQ_rare_factor')
       parameters.RARE=str2num(char(file_acqp(i+1))); %rare factor
   end
   
   if strcmp(file_acqp(i),'##$ACQ_O2_list')     % list of ppm values
       i=i+1;
       %reference freq. here is the water frequency so B0_in_MHz       
       offset = (str2num(char(file_acqp(i+1)))/parameters.B0_in_MHz); 
       parameters.ppm_values = horzcat(parameters.ppm_values,offset);
       point = size(parameters.ppm_values);
       i=i+1;
       while point(2) < parameters.number_of_ppm_values
           %reference freq. here is the water frequency so B0_in_MHz
           offset =(str2num(char(file_acqp(i+1)))/parameters.B0_in_MHz);  
           parameters.ppm_values = horzcat(parameters.ppm_values,offset);
           point = size(parameters.ppm_values);
           i=i+1;
       end                   
   end
end

%sort ppm values

%reading the date and time
%not really needed, but nieeeeeh!
for i=1:size_acqp(1)
   if strncmpi(file_acqp(i),'$$',2)
       line_data=char(file_acqp(i));
       parameters.data_date_and_time=line_data(4:22);
       break
   end
end

clear offset point 
end

% %asks which slice to analyse if number of slices is greater than 1
% %it should be on a different part in the code. Since it is specifically
% %for CEST images.
% if number_of_slices > 1
%     analyzed_slice_index=str2num(char(inputdlg('Which slice to analyze?' ,'',1, {'1'} )));  
% else
%     analyzed_slice_index=1;
% end 
