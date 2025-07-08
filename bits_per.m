function [bitsPerSpike, bitsPerSecond] = bits_per(rateMap, occMap)
% Computes spatial information content
% Inputs:
%   rateMap — matrix of firing rates per spatial bin (e.g., Hz)
%   occMap  — matrix of occupancy probabilities per bin (must sum to 1)
% Outputs:
%   bitsPerSpike  — information per spike (Skaggs et al., 1993)
%   bitsPerSecond — information per second

    % Flatten
    r = rateMap(:);
    p = occMap(:);

    % Remove bins with NaN or zero occupancy
    valid = ~isnan(r) & ~isnan(p) & p > 0;
    r = r(valid);
    p = p(valid);

    R = sum(r .* p);  % Overall mean firing rate (expected r)
    if R == 0
        bitsPerSpike = NaN;
        bitsPerSecond = NaN;
        return
    end

    % Avoid log2(0) by replacing zero rates with eps
    r(r==0) = eps;

    % Compute metrics
    ratio = r / R;
    logTerm = log2(ratio);
    bitsPerSpike  = sum(p .* ratio .* logTerm);
    bitsPerSecond = sum(p .* r .* logTerm);
end
