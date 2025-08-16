function f = plotpartialregression(x, y, z, mask)
%PLOTPARTIALREGRESSION  Partial‐regression scatter of y vs x controlling for z
%   f = plotpartialregression(x,y,z,mask) …
%
% ex:
% plotpartialregression(MI_place_noCSUS_all,MI_CSUS15_all,rates_all, ratemask_all)
% plotpartialregression(bitsper_all(1,:), bitsperCSUS_all(1,:), rates_all, ratemask_all)
% plotpartialregression(bitsper_all(2,:), bitsperCSUS_all(2,:), rates_all, ratemask_all)


% --- 1) define a clean index vector 'temp' so you only keep non‐NaNs:


if nargin < 4
  mask = ones(length(x),1);
end  % std threshold

temp = find(mask == 1)';

% pull out the three column‐vectors
x = x(temp)';        % intertrial bps
y = y(temp)';        % CS+Trace bps
z = z(temp)';        % firing rate

temp = find( ...
    ~isnan(x) & ...
    ~isnan(y) & ...
    ~isnan(z)         );


    % pull out the three column‐vectors
    x = x(temp);        % intertrial bps
    y = y(temp);        % CS+Trace bps
    z = z(temp);        % firing rate



% --- 2) sanity‐check lengths
fprintf('N points = %d\n', numel(x));
assert(~isempty(x),'No data left after masking!');
assert(all(isfinite(x)),'x contains NaNs or Infs');
assert(all(isfinite(y)),'y contains NaNs or Infs');
assert(all(isfinite(z)),'z contains NaNs or Infs');

% --- 3) compute residuals of x and y after regressing out z
X = [ones(numel(z),1), z];
beta_x = X\x;
beta_y = X\y;
r1 = x - X*beta_x;   % intertrial residuals
r2 = y - X*beta_y;   % task residuals

% --- 4) scatter
scatter(r1, r2, 25, 'filled');
hold on; lsline;           % add least‐squares line
%axis square;
xlabel('Intertrial residual bps');
ylabel('CS+Trace residual bps');


% --- 5) annotate with partial correlation
[rp, pp] = partialcorr(x,y,z,'Rows','complete');

title(sprintf('Partial r=%.2f, p=%.3f', rp, pp))
