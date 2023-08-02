function f = findbadTTLs(eventsfile)
%NOT IN USE

ttl =  eventsfile;
ttl(:,1) = ttl(:,1)./1000000;
apart = 1/15;

%finds timestamps with either 0 or 4
ts04 = sort(find(ttl(:,2)==0 | ttl(:,2)==4 ));
sysClock = ttl(ts04,1);
frame_numbers = [1:1:length(sysClock)]';
lost_frames_indices = findlostframes(sysClock);

% Step 6: Recompute the corrected timestamps
count=0
while length(lost_frames_indices)>53 & count<2
  for i=1:length(lost_frames_indices)
    index_to_insert = lost_frames_indices(i) % Adjust the insertion index
    corrected_vector = [sysClock(1:index_to_insert); sysClock(index_to_insert)+apart; sysClock(index_to_insert+1:end)];
    sysClock = corrected_vector;
  end
  count = count+1;
lost_frames_indices = findlostframes(sysClock);
end

length(sysClock);
end



function lost_frames_indices = findlostframes(sysClock)
%Step 1: Estimate the time per frame (eTPF)
eTPF = 1/15;
% Step 3: Calculate the clock differential
clock_differential = sysClock(2:end)-sysClock(1:end-1);
% Step 4: Identify frames with clock differential greater than or equal to eTPF*1.5
lost_frames_indices = find(clock_differential>((eTPF)*1.5));
% Step 5: Determine the number of dropped frames
num_dropped_frames = length(lost_frames_indices)
end





%corrected_timestamps = sysClock;
%for i = 1:length(num_dropped_frames)
%    % Correct timestamps after each lost frame
%    current_lost_frame_index=(lost_frames_indices(i))
%    correction = sysClock(current_lost_frame_index):eTPF:sysClock(current_lost_frame_index+1);
%    corrected_timestamps(lost_frames_indices(i)+1:end) = corrected_timestamps(lost_frames_indices(i)+1:end) + correction;
%end



%{
for checking
%finds where interval between TTL 0s is greater than expected
ts0 = (find(ttl(:,2)==0));
diff = ttl(ts0(2:end,1))-ttl(ts0(1:end-1,1));
missing0 = (find(diff>(apart*2*1.5)));
missing0 = missing0*2

%finds where interval between TTL 4s is greater than expected
ts4 = (find(ttl(:,2)==4));
diff = ttl(ts4(2:end,1))-ttl(ts4(1:end-1,1));
missing4 = (find(diff>(apart*2*1.5)));
missing4 = missing4*2
%}
