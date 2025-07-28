function [bitsPerSpike, bitsPerSecond] = bits_perCSUS(rateMap, occMap)
% Computes spatial information content
%bits/spike: Sum of (occprobs * mean firing rate per bin / meanrate) * log2 (mean firing rate per bin / meanrate)
%bits/sec: Sum of (occprobs * mean firing rate per bin) * log2 (mean firing rate per bin / meanrate)
% Inputs:
%   rateMap — matrix of firing rates per spatial bin (e.g., Hz)
%   occMap  — matrix of occupancy probabilities per bin (must sum to 1)
% Outputs:
%   bitsPerSpike  — information per spike (Skaggs et al., 1993)
%   bitsPerSecond — information per second


% -------------------------------------------------------------------------
% 1. sanity checks
% -------------------------------------------------------------------------
if isequal(size(rateMap), size(occMap))==0
  rateMap = rateMap';
end
if isequal(size(rateMap), size(occMap))==0
  size(rateMap)
  size(occMap)
end
assert(isequal(size(rateMap), size(occMap)), ...
    'rateMap and occMap must be the same size');
occSum = nansum(occMap(:));
assert(abs(occSum - 1) < 1e-6, 'occMap must sum to 1 (got %.3g)', occSum);



pOcc      = occMap;     % column vector
rBin      = rateMap;        % column vector



% -------------------------------------------------------------------------
% 3. compute information terms (vectorised)
% -------------------------------------------------------------------------
% 3. compute information terms (vectorised, skip zero-rate bins)
meanRate = nanmean(rBin);                 % overall mean firing (Hz)

posMask  = rBin > 0;                      % ignore bins with r=0
pOccPos  = pOcc(posMask);
rPos     = rBin(posMask);
ratio    = rPos ./ meanRate;              % r_i / r̄   for positive-rate bins

bitsPerSpike  = sum( pOccPos .* ratio .* log2(ratio) );
bitsPerSecond = sum( pOccPos .* rPos   .* log2(ratio) );

end
