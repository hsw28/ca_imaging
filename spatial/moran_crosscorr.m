function [all_info forward_info backward_info] = moran_crosscorr(currentcellcenter, currentfieldcenter, spike_times, pos)


set(0,'DefaultFigureVisible', 'off');

spike_times = spike_times';
currentfieldcenter = currentfieldcenter';


  %finds fast times
  velthreshold = 12;
  vel = ca_velocity(pos);
  vel(1,:) = smoothdata(vel(1,:), 'gaussian', 30.0005); %originally had this at 30, trying with 15 now
  badvel = find(vel(1,:)>=velthreshold); %INDICES OF BAD VELS
  goodpos = pos(:,1:2);
  goodpos(badvel, 2) = NaN;
  %gets direction of fast times
  fwd = [];
  bwd = [];
  for z = 1:length(goodpos)
    if isnan(goodpos(z,2))==0
          if goodpos(max(z-15, 1),2)-goodpos(min(z+15,length(pos)),2)>0
              fwd(end+1) = (z); %indexes for forward
          else
              bwd(end+1) = (z); %indexes for backwards
          end
      end
    end



  goodpos(fwd, 2) = 1;
  goodpos(bwd, 2) = 2; %now have pos with times and then a NaN for not moving, 1 for moving east, 2 for moving west



  fwd_index = find(goodpos(:, 2)==1);
  bwd_index = find(goodpos(:, 2)==2);

  fwd_time = goodpos;
  fwd_time(:,2) = NaN;
  fwd_time(fwd_index,2) = 1;

  bwd_time = goodpos;
  bwd_time(:,2) = NaN;
  bwd_time(bwd_index,2) = 2;

  n = 1;
  startindex = [];
  endindex = [];
  while n<=length(fwd_time)
    if n<=length(fwd_time) & (fwd_time(n,2))==1
      startindex(end+1) = n;
      while n<length(fwd_time) & (fwd_time(n,2))==1
        n = n+1;
      end
    endindex(end+1) = n;
  end
    n = n+1;
  end
  fwd_ind = [startindex; endindex]; %index of forwards
  length_fwd = endindex-startindex;
  good_fwd = find(length_fwd > 30); %find forward times over 1 second long
  fwd_ind = fwd_ind(:,good_fwd);
  % (fwd_ind) is [2, n] where n is thee start and stop times

  n = 1;
  startindex = [];
  endindex = [];

  while n<=length(bwd_time)
    if n<=length(bwd_time) & (bwd_time(n,2))==2
      startindex(end+1) = n;
      while n<length(bwd_time) & (bwd_time(n,2))==2
        n = n+1;
      end
    endindex(end+1) = n;
  end
    n = n+1;
  end
  bwd_ind = [startindex; endindex]; %index of forwards
  length_bwd = endindex-startindex;
  good_bwd = find(length_bwd > 30); %find forward times over 1 second long
  bwd_ind = bwd_ind(:,good_bwd);
  % (fwd_ind) is [2, n] where n is the index of start and stop times


  pos = pos';




  vectorsize = length(currentcellcenter).*(length(currentcellcenter)-1);
  %forward
  count = 1;
  dir = NaN(1, vectorsize, 'single'); %direction, 1 for east/for, 2 for west/back
  index1 = NaN(1, vectorsize, 'single'); %cell index 1
  index2 = NaN(1, vectorsize, 'single'); %cell index 2
  field1 = NaN(1, vectorsize, 'single'); %field center cell 1
  field2 = NaN(1, vectorsize, 'single'); %field center cell 2
  field_dist = NaN(1, vectorsize, 'single'); %difference in field centers
  cell_dist = NaN(1, vectorsize, 'single'); %distance bbetween cell centers ((mm))
  max_cc = NaN(1, vectorsize, 'single'); %max cross correelation
  time_cc = NaN(1, vectorsize, 'single'); %offset
  corr_max_av = NaN(1, vectorsize, 'single');
  corr_time_av = NaN(1, vectorsize, 'single');
  lineup = NaN(1, vectorsize, 'single');
for k=1:length(currentcellcenter)-1
    for j = 2:length(currentcellcenter)
      if count<0
      numlineup = 0;

    center_fwd_1 = currentfieldcenter(1,k).*4;
    center_fwd_2 = currentfieldcenter(1,j).*4;

    cell1_center = currentcellcenter(1:2, (k));
    cell2_center = currentcellcenter(1:2, (j));
    dis = abs(norm(cell1_center-cell2_center)).*0.0055; %distance between points in pixels,convert to mm

    cell1_spikes = spike_times(:,k);
    cell2_spikes = spike_times(:,j);

    corr_max = NaN(length(fwd_ind),1);
    corr_time = NaN(length(fwd_ind),1);

      for n=1:length(fwd_ind)

        start_time = pos(1, fwd_ind(1,n));
        end_time = pos(1, fwd_ind(2,n));

        cell1_spikes_in = intersect(find(cell1_spikes<end_time), find(cell1_spikes>start_time));
        cell2_spikes_in = intersect(find(cell2_spikes<end_time), find(cell2_spikes>start_time));
        cell1_spikes_in = cell1_spikes(cell1_spikes_in);
        cell2_spikes_in = cell2_spikes(cell2_spikes_in);
        if length(cell1_spikes_in)>3 & length(cell2_spikes_in)>3
          cell1_spikes_in = histcounts(cell1_spikes_in, start_time:.01:end_time);
          cell2_spikes_in = histcounts(cell2_spikes_in, start_time:.01:end_time);
          corrs = Hcorr(cell1_spikes_in, cell2_spikes_in);
          [val, inx] = max(corrs(1,:));
          corr_max(n) = val;
          corr_time(n) = (corrs(2,inx)).*.01;
          numlineup = numlineup+1;
        end
      end



    dir(count) = 1; %direction, 1 for east/for, 2 for west/back
    index1(count)= k; %cell index 1
    index2(count) = j; %cell index 2
    field1(count) = center_fwd_1; %field center cell 1
    field2(count) = center_fwd_2; %field center cell 2
    field_dist(count) = abs(center_fwd_1-center_fwd_2); %diifference in field centers
    cell_dist(count) = abs(norm(cell1_center-cell2_center)); %distance bbetween cell centers ((mm))
    corr_max_av(count) = nanmean(corr_max); %max cross correelation
    corr_time_av(count) = nanmean(corr_time); %offset
    lineup(count) = numlineup;
    count = count+1;
    if rem(count,1000)==0
      count;
    end

  else
    count = count+1;
  end

  end
end
count
forward_info = [dir; index1; index2; field1; field2; field_dist; cell_dist; corr_max_av; corr_time_av; lineup];


%back
count = 1;
dir = NaN(1, vectorsize, 'single'); %direction, 1 for east/for, 2 for west/back
index1 = NaN(1, vectorsize, 'single'); %cell index 1
index2 = NaN(1, vectorsize, 'single'); %cell index 2
field1 = NaN(1, vectorsize, 'single'); %field center cell 1
field2 = NaN(1, vectorsize, 'single'); %field center cell 2
field_dist = NaN(1, vectorsize, 'single'); %difference in field centers
cell_dist = NaN(1, vectorsize, 'single'); %distance bbetween cell centers ((mm))
max_cc = NaN(1, vectorsize, 'single'); %max cross correelation
time_cc = NaN(1, vectorsize, 'single'); %offset
corr_max_av = NaN(1, vectorsize, 'single');
corr_time_av = NaN(1, vectorsize, 'single');
lineup = NaN(1, vectorsize, 'single');
for k=1:length(currentcellcenter)-1
  for j = 2:length(currentcellcenter)
    if count>=50000 & count<50000000


  numlineup = 0;
  center_fwd_1 = currentfieldcenter(2,k).*4;
  center_fwd_2 = currentfieldcenter(2,j).*4;

  cell1_center = currentcellcenter(1:2, (k));
  cell2_center = currentcellcenter(1:2, (j));
  dis = abs(norm(cell1_center-cell2_center)).*0.0055; %distance between points in pixels

  cell1_spikes = spike_times(:,k);
  cell2_spikes = spike_times(:,j);

  corr_max = NaN(length(bwd_ind),1);
  corr_time = NaN(length(bwd_ind),1);

    for n1=1:length(bwd_ind)

      start_time = pos(1, bwd_ind(1,n1));
      end_time = pos(1, bwd_ind(2,n1));

      cell1_spikes_in = intersect(find(cell1_spikes<end_time), find(cell1_spikes>start_time));
      cell2_spikes_in = intersect(find(cell2_spikes<end_time), find(cell2_spikes>start_time));
      cell1_spikes_in = cell1_spikes(cell1_spikes_in);
      cell2_spikes_in = cell2_spikes(cell2_spikes_in);
      if length(cell1_spikes_in)>3 & length(cell2_spikes_in)>3
        cell1_spikes_in = histcounts(cell1_spikes_in, start_time:.01:end_time);
        cell2_spikes_in = histcounts(cell2_spikes_in, start_time:.01:end_time);
        corrs = Hcorr(cell1_spikes_in, cell2_spikes_in);
        [val, inx] = max(corrs(1,:));
        corr_max(n1) = val;
        corr_time(n1) = (corrs(2,inx))*.01;
        numlineup = numlineup+1;
      end
    end



  dir(count) = 2; %direction, 1 for east/for, 2 for west/back
  index1(count)= k; %cell index 1
  index2(count) = j; %cell index 2
  field1(count) = center_fwd_1; %field center cell 1
  field2(count) = center_fwd_2; %field center cell 2
  field_dist(count) = abs(center_fwd_1-center_fwd_2); %diifference in field centers
  cell_dist(count) = abs(norm(cell1_center-cell2_center)); %distance bbetween cell centers ((mm))
  corr_max_av(count) = nanmean(corr_max); %max cross correelation
  corr_time_av(count) = nanmean(corr_time); %offset
  lineup(count) = numlineup;
  count = count+1;
  if rem(count,5000)==0
    count
  end

else
  count = count+1;
end




end
end
count

backward_info = [dir; index1; index2; field1; field2; field_dist; cell_dist; corr_max_av; corr_time_av; lineup];

all_info = [forward_info, backward_info];

%try seperating
%put in markers to see where quits
