%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Afroditi Eleftheriou                                                    %
% 25.11.2017                                                              %
%                                                                         %
% Checks if the bruker numbering of scans is continuous or not;           %
% If numbering is not continuous, user has to type the exact order.       %
% Input : first                                                           %
%         number_of_scans                                                 %
%         regular_order                                                   %
% Output: scan (that has the bruker numbers of scans in order)            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[scan] = scan_ordering(first,number_of_scans,regular_order)

% checks if the numbering of scans is continuous or not
if regular_order == 1 | strcmp(regular_order, 'yes')| strcmp(regular_order,'Yes')| ...
   strcmp(regular_order, 'YES')| strcmp(regular_order, 'Y')| strcmp(regular_order, 'y')
   
      scan(1) = first;
      for i=2:number_of_scans
          scan(i) = scan(i-1)+1; %loads 'scan' with the continuous numbering
      end
else
      prompt = {'Put the numbers of the scans folders in order, separated by space (e.g. 10 12 13 17)'};
      defaults = { '10 12 13 17'}; 
      dlg_title = 'Better if ordered... Don''t scew up again!';
      num_lines = 1;
      scan_numbers = inputdlg(prompt,dlg_title,num_lines,defaults); %takes user input individual numbers of scans 
      scan_numbers = str2num(scan_numbers{:});
      scan = scan_numbers; %loads 'scan' with the user defined numbering
end
end