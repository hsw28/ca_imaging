function f  = interp_pos(pos)

if size(pos,1)<size(pos,2)
pos = pos';
end

X = pos(:,1);
X = inpaint_nans(X, 4);

Y = pos(:,2);
Y = inpaint_nans(Y, 4);

f = [X,Y];
