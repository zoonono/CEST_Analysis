i_time = 1;
flag = 0;
for ii = 1 : input.number_of_scans
  ppm_number = 0;  
  for jj = 1 : cest_scan(ii).number_of_ppm_values
    if cest_scan(ii).ppm_values(jj) == cest_scan(ii).ppm_values(1) && jj ~= 1    
      ppm_number = 1; 
      timeline(i_time,ppm_number) = cest_scan(ii).ppm_values(jj);
      data_in_ppm_sets(i_time) = cest_scan(ii);
      i_time = i_time + 1;
    else
      ppm_number = ppm_number + 1;
      timeline(i_time,ppm_number) = cest_scan(ii).ppm_values(jj);
      data_in_ppm_sets(i_time) = cest_scan(ii);
    end
  end  
  i_time = i_time + 1;
end

for i = 1 : 270
  data_in_ppm_sets(i).timeline(:) = timeline(i,:);
end
