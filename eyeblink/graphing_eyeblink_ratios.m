function f = graphing_eyeblink_ratios(structure)
  %makes a matrix of averaged response values for a full day

  fields = fieldnames(structure);

  values = NaN(3,length(fields));

  for i = 1:numel(fields)
    fieldName = fields{i};
    current = structure.(fieldName);


    values(:,i)= (sum(current)./sum(current(:,1)));
  end

f = values;
