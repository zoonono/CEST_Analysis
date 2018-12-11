%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Afroditi Eleftheriou                                                    %
% 14.12.2017                                                              %
%                                                                         %
% Print data in txt file, with headers                                    %
%                                                                         %
%  Input: ppm_values                                                      %
%         number_of_scans                                                 %
%         scan_data(#scans,#ppms)                                         %
%         filename                                                        %
%         path                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[text_file] = print_in_txt(ppm_values,number_of_scans, scan_data,filename,path)

header = string(zeros(1,number_of_scans));
%create header: scan1 scan2 etc...
for i = 1:number_of_scans
  header(i) = strcat('scan',num2str(i));
end

%w+: Open or create new file for reading and writing. Discard existing contents, if any.
fid = fopen([path filename], 'w+');

fprintf(fid,'%7s','ppm');
fprintf(fid,'%15s\t', header{:});
fprintf(fid,'\n');

for j = 1:length(ppm_values)
  fprintf(fid, '%7.2f\t', ppm_values(j));
  %scan_data(#scan,#ppm)
  for i = 1:number_of_scans
    fprintf(fid, '%15.5f\t', scan_data(i,j));
  end
  fprintf(fid,'\n');
end
fclose(fid);

text_file = filename;

end