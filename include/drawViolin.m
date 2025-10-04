function drawViolin(x0, data, width)
  if nargin<3 || isempty(width), width = 0.28; end   % <- narrower default
  data = data(~isnan(data));
  if numel(data)<3
      plot(x0+0.02*randn(size(data)), data, 'k.'); return;
  end
  [f,y] = ksdensity_fallback(data);
  f = f / max(f) * width;  % <- use width
  patch([x0 - f, fliplr(x0 + f)], [y, fliplr(y)], [0.6 0.7 0.9], ...
        'EdgeColor',[0.3 0.3 0.5], 'FaceAlpha',0.6, 'LineWidth',1);
end
