%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Afroditi Eleftheriou                                                    %
% 23.11.2017                                                              %
%                                                                         %
% Reads 'method' file                                                     %
%                                                                         %
% Input: folder_number                                                    %
%        directory                                                        %
%                                                                         %
% Output: parameters.B1_power                                             %
%         parameters.saturation_pulse_duration                            %
%         parameters.NAverages                                            %
%         parameters.orientation_read_direction                           %
%         parameters.orientation_slice                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These two commented lines are only for testing
%
% folder_number = 17; 
% directory = 'D:\CEST_data\MRI data\mouse_3OMGCEST\20171031_A_3OMGCEST';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[parameters] = read_method(folder_number,directory)

%sets the directory of the scan folder (the name of the folder is a number)
folder_number = num2str(folder_number); %the index that is being looped for all scans in the main program
directory_scan = strcat(directory,'\',folder_number); %directory of the scan
cd(directory_scan);

file_method = textread('method','%s','delimiter','=','whitespace','');
size_method = size(file_method);

for i=1:size_method(1)
    if strcmp(file_method(i),'##$PVM_MagTransPower')
        parameters.B1_power = str2num(char(file_method(i+1))); % B1 power in microTesla
    end
    if strcmp(file_method(i),'##$PVM_MagTransModuleTime')
       parameters.saturation_pulse_duration = str2num(char(file_method(i+1))); %saturation impulse duration in sec
    end
    if strcmp(file_method(i),'##$PVM_NAverages')
       parameters.NAverages = str2num(char(file_method(i+1))); % Number of Averages
    end    
    if strcmp(file_method(i),'##$PVM_SPackArrReadOrient')
        % finds the image readout orientation: normal L_R    
        parameters.orientation_readout = (char(file_method(i+2))); %orientation of read diraction
    end
    if strcmp(file_method(i),'##$PVM_SPackArrSliceOrient')
        parameters.orientation_slice = (char(file_method(i+2))); %orientation of dlice diraction
    end
end
end
