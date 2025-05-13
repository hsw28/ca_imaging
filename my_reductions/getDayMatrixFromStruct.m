function [X, y] = getDayMatrixFromStruct(animal, dateStr, win, nBins, Fs)
    csVar = ['CS_' dateStr];
    caVar = ['CA_traces_' dateStr];
    ftVar = ['CA_time_' dateStr];
    ttVar = ['EMGts_' dateStr];

    if ~isfield(animal.CS_times, csVar) || ~isfield(animal.Ca_traces, caVar) || ~isfield(animal.Ca_ts, ftVar) || ~isfield(animal.CRs, ttVar)
        X = []; y = [];
        return
    end

    CSon = animal.CS_times.(csVar);
    F = animal.Ca_traces.(caVar);
    t = animal.Ca_ts.(ftVar);
    t = t(:,2)/1000; %changes from miliseconds
    t = t(1:2:end); % samples every other timepoint bc Ca_ts is 15hz and F is 7.5hz
    t = t(1:size(F,2)); %sometimes sampling every other timepoint makes Ca_ts one longer than F, so truncate to time of F
    y = animal.CRs.(ttVar);

    %fprintf('Raw y labels for %s: %d zeros, %d ones, %d NaNs\n', dateStr, sum(y==0), sum(y==1), sum(isnan(y)));

    if isempty(CSon) || numel(CSon) < 2 || all(isnan(CSon))
        X = []; y = [];
        return
    end

    nC = size(F,1);
    nT = numel(CSon);
    X = nan(nC, nBins, nT, 'single');
    usedTrialIdx = false(nT,1);

    for k = 1:nT
        t0 = CSon(k) + win(1);
        t1 = CSon(k) + win(2);

        if t0 < t(1) || t1 > t(end)
            fprintf('Skipping trial %d: outside valid time range (t0 = %.2f, t1 = %.2f)\n', ...
                k, t0, t1);
            continue;
        end

        idx = find(t >= t0, 1, 'first') + (0:nBins-1);
        if any(idx > length(t))
            fprintf('Skipping trial %d: would exceed trace length.\n', k);
            continue;
        end
        trace = F(:, idx);
        

        trace = F(:, idx);

        if size(trace,2) < nBins
            fprintf('⚠️ Padding trial %d (CSon = %.2f): only %d frames, expected %d\n', ...
                k, CSon(k), size(trace,2), nBins);
            trace(:, end+1:nBins) = NaN;  % <-- This is what's giving you NaNs
        elseif size(trace,2) > nBins
            trace = trace(:, 1:nBins);  % optional trimming
        end

        X(:,:,k) = trace;
        usedTrialIdx(k) = true;
    end


    if ~isempty(X)
        reshaped = reshape(X, nC, []);
        mu = mean(reshaped, 2, 'omitnan');
        X = bsxfun(@minus, X, mu);
    end
end
