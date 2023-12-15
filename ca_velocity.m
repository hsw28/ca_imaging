function v = ca_velocity(pos);
%computes velocity from position data as velocity(pos, varargin). if you put a 1 in varargin, velocity will be directional from junction of forced arms (380, 360)
%computes velocity. input a [#ofpoints, 3] vector, where first column is time, second is x, third is y
%linear track is 248.92cm and about 620 pixels = 2.5pixels per cm

file = pos';
%file = fixpos(pos);


t = file(1, :);
xpos = (file(2, :))';
ypos = (file(3, :))';



velvector = [];
timevector = [];

s = size(t,2);

for i = 2:s-1
	%find distance travelled
	if t(i)~=t(i-1)
		hypo = hypot((xpos(i-1)-xpos(i+1)), (ypos(i-1)-ypos(i+1)));
		vel = hypo./((t(i+1)-t(i-1)));
		velvector(end+1) = vel;
		timevector(end+1) = t(i);
	end
end

%velvector = filloutliers(velvector, 'pchip', 'movmedian',10);

%v = hampel(velvector, 30, 3);
v = smoothdata(velvector,'gaussian',7);
v = v(1:length(timevector));

v = [(v/1.000); timevector];
