%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Afroditi Eleftheriou                                                    %
% 25.11.2017                                                              %
%                                                                         %
% Reads 'reco' file                                                       %
%                                                                         % 
% Input: folder_number                                                    %
%        directory                                                        %
%                                                                         %
% Output: parameters.reco_bit      (uint16 or uint32)                     %
%         parameters.size_x        in mm                                  %
%         parameters.size_y        in mm                                  %
%         parameters.slope_values  for each ppm (Repetition)              %
%         parameters.slope         only one (they're the same anyway)     %
%         parameters.FOV           in mm                                  %
%         parameters.reco_wordtype                                        %
%         parameters.directory_reco_and_2dseq                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These two commented lines are only for testing
%
% folder_number = 17; 
% directory = 'D:\CEST_data\MRI data\mouse_3OMGCEST\20171031_A_3OMGCEST';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[parameters] = read_reco(folder_number,directory)

% goes to the directory of the scan folder (the name of the folder is a number)
folder_number = num2str(folder_number); %folder_number is the index that is being looped for all scans in the main program
directory_scan = strcat(directory,'\',folder_number); %directory of the scan
directory_reco = strcat(directory_scan,'\pdata\1');
cd(directory_reco)

parameters.directory_reco_and_2dseq = directory_reco;

file_reco=textread('reco','%s','delimiter','=','whitespace','');
size_reco=size(file_reco);

parameters.slope_values = [];

for i=1:size_reco(1)
    if strcmp(file_reco(i),'##$RECO_fov')
        %flag_FOV = i; %comment out, only for test
        parameters.FOV = str2num(char(file_reco(i+2))); % in cm
        parameters.FOV = parameters.FOV * 10;      % in mm
    end    
    if strcmp(file_reco(i),'##$RECO_size')
        %flag_size_matrix=i; %comment out, only for test
        size_matrix = str2num(char(file_reco(i+2)));
        parameters.size_x = size_matrix(2);  %not sure why it is second and not first
        parameters.size_y = size_matrix(1); 
    end    
    if strcmp(file_reco(i),'##$RECO_wordtype')
        %flag_wordtype = i; %comment out, only for test
        parameters.reco_wordtype = (char(file_reco(i+1)));
    end
    if strcmp(file_reco(i),'##$RECO_map_slope')
        %flag_slope = i; %comment out, only for test
        number_of_slope_values = str2num(char(file_reco(i+1)));
        
        %extracts all #slope_index slopes. 
        %There are, most probably, all the same.
        i=i+1;
        slopes_in_one_line=(str2num(char(file_reco(i+1))));                    
        parameters.slope_values=horzcat(parameters.slope_values,slopes_in_one_line);
        size_of_slope_values = size(parameters.slope_values);
        i=i+1;
        while size_of_slope_values(2) < number_of_slope_values
            slopes_in_one_line =(str2num(char(file_reco(i+1))));  
            parameters.slope_values = horzcat(parameters.slope_values,slopes_in_one_line);
            size_of_slope_values = size(parameters.slope_values);
            i=i+1;
        end        
    end
end

if strcmp(parameters.reco_wordtype,'_32BIT_SGN_INT')
    parameters.reco_bit='uint32';
elseif strcmp(parameters.reco_wordtype,'_16BIT_SGN_INT')
    parameters.reco_bit='uint16';
end

parameters.slope = parameters.slope_values(1);

clear number_of_slope_values size_of_slope_values slopes_in_one_line
clear directory_reco directory_scan
end