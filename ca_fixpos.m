function f = ca_fixpos(pos)
%converts position file to time, x, y and interpolates missing values
%REMEMBER TO GET PIXEL TO INCH CONVERSION

length_in_seconds = length(pos)./30;

pos = pos';
time = (0:(length_in_seconds./length(pos)):length_in_seconds);
%x = inpaint_nans(pos(1,:), 2);
%y = inpaint_nans(pos(2,:), 2);

y = inpaint_nans(pos, 2);

%pos = [time(2:end); x; y];
pos = [time(2:end); y];


f = pos';
