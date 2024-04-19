function f = DEPca_fixpos(pos, timestamps)
%converts position file to time, x, y and interpolates missing values
%USE dlc_fixpos for deeeplabcut videoos


%DELETE NEXT LINES
%startpos = pos(1,1)*30;
%timestamps = timestamps(startpos:end);
%pos = pos(:,[2,3]);


pos = pos';
time = timestamps/1000;
%x = inpaint_nans(pos(1,:), 2);
%y = inpaint_nans(pos(2,:), 2);


for k=1:length(pos)
  if isnan(pos(2,k))==0
    pos=pos(:,k:end);
    time = time(k:end);
    break
  end
end


y = inpaint_nans(pos, 2);

%pos = [time(2:end); x; y];

time = time(1:length(y));
pos = [time(1:end)'; y];


vel = ca_velocity(pos');

vel = smoothdata(vel, 'gaussian', 30);

pos = pos';


highvstart = find(vel(1,1:5400)>100);


stop = length(pos)-2700;
highvend = find(vel(1,stop:end)>100);
pos = pos(highvstart:max(highvend)+stop,:);


f = pos;
