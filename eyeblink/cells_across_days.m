function alignment_response = cells_across_days(ratios_structure, alignment)

  %for every cell, finds the pre-trial response, CS response, and US response across all days the cell appears
  %ratios is from eyeblink_ratios(times_US, Ca_peaks), where ratios are [pretrial_sum, CS_sum, US_sum, CS_change, US_change];

fields = fieldnames(ratios_structure);
pre_response = NaN(size(alignment,1),size(alignment,2));
CS_resp = NaN(size(alignment,1),size(alignment,2));
US_resp = NaN(size(alignment,1),size(alignment,2));
CS_ratio = NaN(size(alignment,1),size(alignment,2));
US_ratio = NaN(size(alignment,1),size(alignment,2));


for i = 1:size(alignment,1)
  for j = 1:size(alignment,2)

      %day we want is index j
      %cell we want is whatever is in i
      fieldName = fields{j};
      current_day = ratios_structure.(fieldName);
      cell_number = alignment(i,j);
      if cell_number>0 & length(current_day)>1
          pre_response(i,j) = current_day(cell_number,1);
          CS_resp(i,j) = current_day(cell_number,2);
          US_resp(i,j) = current_day(cell_number,3);
          CS_ratio(i,j)= current_day(cell_number,4);
          US_ratio(i,j)= current_day(cell_number,5);

      end
  end
end

alignment_response = struct();
alignment_response.pre_response = pre_response;
alignment_response.CS_resp = CS_resp;
alignment_response.US_resp = US_resp;
alignment_response.CS_ratio = CS_ratio;
alignment_response.US_ratio = US_ratio;
