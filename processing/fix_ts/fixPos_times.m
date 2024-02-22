
function fixed_pos = fixPos_times(pos)

fields = fieldnames(pos);
for i = 1:numel(fields)
    fieldName = fields{i};
    current_pos = pos.(fieldName);

  if size(current_pos,2)>size(current_pos,1)
    current_pos = current_pos';
  end



  if contains(fieldName, 'oval')
    for_rec_1_for_oval_2 = 2;
  else
    for_rec_1_for_oval_2 = 1;
  end


  if current_pos(1,1)>1 || current_pos(2,1)>1
    current_pos(:,1) = current_pos(:,1)/ 1000;
  end

  xpos = current_pos(:,2);
  ypos = current_pos(:,3);



  if for_rec_1_for_oval_2 == 2 %oval
    xpos = xpos*.14;
    xpos = xpos-min(xpos);
    ypos = ypos*.15;
    ypos = ypos-min(ypos);
  else
      xpos = xpos*.19;
      xpos = xpos-min(xpos);
      ypos = ypos*.15;
      ypos = ypos-min(ypos);
  end


  current_pos(:,2) = xpos;
  current_pos(:,3) = ypos;

pos.(fieldName) = current_pos;
end

fixed_pos = pos;



fprintf('positions fixed')
