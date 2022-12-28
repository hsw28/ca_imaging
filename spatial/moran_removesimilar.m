function [f count] = moran_removesimilar(cellcenter, fieldcenter, distance_in_mm)
% removes any cells that might be duplicates based on their distance apart (distance_in_mm)
% and similar place fields in both directions

similar_fields = 10; %distance apart in cm fields need to be less than in order to be similar

set(0,'DefaultFigureVisible', 'off');
if isstruct(cellcenter)==1
  cellnames = fieldnames(cellcenter);
  fieldcenternames = fieldnames(fieldcenter);
  daynum = length(cellnames);
else
  daynum =1;
end

count = [];

for z=1:(daynum)
  if isstruct(cellcenter)==1
    name = char(cellnames(z));
    currentcellcenter = cellcenter.(name);
    name = char(fieldcenternames(z));
    currentfieldcenter = fieldcenter.(name);
  else
    currentcellcenter = cellcenter;
    currentfieldcenter = fieldcenter;
  end


  center_fwd = currentfieldcenter(:,1).*4;
  center_bwd = currentfieldcenter(:,2).*4;

  %differeences in field centers
  diff_fwd = abs(bsxfun(@minus, center_fwd, center_fwd'));
  diff_bwd = abs(bsxfun(@minus, center_bwd, center_bwd'));

  ind_fwd = find(diff_fwd(:) < similar_fields);
  ind_bwd = find(diff_bwd(:) < similar_fields);

  potential_same = intersect(ind_fwd, ind_bwd);


  [cell1_index,cell2_index] = ind2sub([length(center_fwd), length(center_fwd)],potential_same);

  num = [];
  same1 = [];
  same2 = [];
  for k=1:length(cell1_index)
    cell1_center = currentcellcenter(1:2, cell1_index(k));
    cell2_center = currentcellcenter(1:2, cell2_index(k));
    dis = abs(norm(cell1_center-cell2_center)); %distance between points in pixels
    dis = dis*0.0055; %convert to mm
    if dis<distance_in_mm & (cell1_index(k)~=cell2_index(k))
      %a = (find(same2==cell1_index(k)));
      %b = (find(same1==cell2_index(k)));
      %if length(a)~=length(b) | length(a-b)==0 | (a-b)~=0
      %zz = (ismember(num,[cell2_index(k); cell1_index(k)]))
      if sum(ismember(num,[cell2_index(k); cell1_index(k)]))<1
        same1(end+1) = cell1_index(k);
        same2(end+1) = cell2_index(k);
        num = [same1; same2];
      end
    end
  end

    num = [same1; same2]
    nouse = randi(1:2,length(same1),1);
    nowant = [];
    for z=1:length(nouse)
      rep1 = length(find(num(:)==num(1,z)));
      rep2 = length(find(num(:)==num(2,z)));
      if rep1<2 & rep2<2
        nowant(end+1) = num(nouse(z), z);
      end
      if rep1>1
        nowant(end+1) = num(1, z);
      end
      if rep2>2
        nowant(end+1) = num(2, z);
      end
    end

    if isstruct(cellcenter)==1
    new_center.(name) = currentcellcenter(:, setdiff(1:end,nowant));
    new_field.(name) = currentfieldcenter(setdiff(1:end,nowant), :);
    else
      new_center = currentcellcenter(:, setdiff(1:end,nowant));
      new_field = currentfieldcenter(setdiff(1:end,nowant), :);
    end

    count(end+1) = length(nowant);


end


percents = shuffled_moran(new_center, new_field, 1000);
f = percents;
