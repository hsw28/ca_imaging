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
  %    fprintf('oval')
  %  end



    % Calculate the center of the traveled area
    minx = min(current_pos(:,2));
    miny = min(current_pos(:,3));
    maxx = max(current_pos(:,2));
    maxy = max(current_pos(:,3));
    center_x = (maxx+minx)/2;
    center_y = (maxy+miny)/2;

    % Calculate the shift required to move the center to (a, b)
    shift_x = a - center_x;
    shift_y = b - center_y;

    % Apply the shift to all positions
    shifted_positions = [current_pos(:, 1), current_pos(:, [2, 3]) + repmat([shift_x, shift_y], size(current_pos, 1), 1)];


pos.(fieldName) = shifted_positions;

end

f = pos;
