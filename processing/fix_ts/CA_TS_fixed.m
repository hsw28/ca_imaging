function CA_TS_fixed = fixCA_times(CA_TS)
fields = fieldnames(CA_TS);
  for i = 1:numel(fields)
    fieldName = fields{i};
    current_CA = CA_TS.(fieldName);
  if current_CA(1)>1 || current_CA(2)>1
    current_CA(1,:) = current_CA(1,:)/ 1000;
  end
CA_TS.(fieldName) = current_CA;
end
CA_TS_fixed = CA_TS;
fprintf('calcium times fixed')
end
