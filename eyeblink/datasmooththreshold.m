function f = datasmooththreshold(Ca_trace)

%first denoised ΔF/F by smoothing with a Gaussian filter with a length of 5 bins
%thresholded the result so that values less than 2 robust standard deviations across the time series were set to 0.

outtrace = NaN(size(Ca_trace));
for k = 1:size(Ca_trace, 1)
  curtrace = Ca_trace(k,:);
  w = gausswin(5);
  smoothed = filter(w,1,curtrace);
  dev = 2*mad(smoothed);
  av = nanmean(smoothed);
  setzero = find(smoothed<(av-dev));
  smoothed(setzero) = 0;
  outtrace(k,:) = smoothed;
end


f = outtrace;
