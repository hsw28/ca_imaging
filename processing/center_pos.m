function f = center_pos(pos)


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

    a = 100; %center
    b = 100; %center

    %%comment me out
  %  if for_rec_1_for_oval_2 == 1
  %    current_pos(:,2) = (current_pos(:,2)).*1.1538;
  %  end
  %  if for_rec_1_for_oval_2 == 2
  %    current_pos(:,3) = (current_pos(:,3)).*1.06;
  %    current_pos(:,3) = (current_pos(:,3)).*1.13;
  %  end



    % Calculate the center of the traveled area
    center_x = nanmedian(current_pos(:,2));
    center_y = nanmedian(current_pos(:,3));

    % Calculate the shift required to move the center to (a, b)
    shift_x = a - center_x;
    shift_y = b - center_y;

    % Apply the shift to all positions
    shifted_positions = [current_pos(:, 1), current_pos(:, [2, 3]) + repmat([shift_x, shift_y], size(current_pos, 1), 1)];


pos.(fieldName) = shifted_positions;

end

f = pos;
