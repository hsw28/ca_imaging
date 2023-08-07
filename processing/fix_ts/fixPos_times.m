
function fixed_pos = fixPos_times(pos)
fields = fieldnames(pos);
for i = 1:numel(fields)
    fieldName = fields{i};
    current_pos = pos.(fieldName);
  if current_pos(1,1)>1 || current_pos(2,1)>1
    current_pos(:,1) = current_pos(:,1)/ 1000;
  end
pos.(fieldName) = current_pos;
end
fixed_pos = pos;
fprintf('positions fixed')
end
