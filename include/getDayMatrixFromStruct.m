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
    t = animal.Ca_ts.(ftVar); t = t(:,2)/1000; t = t(1:2:end); t = t(1:size(F,2));
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
        if t0 < t(1) || t1 > t(end), continue; end

        idx = find(t >= t0 & t < t1);
        if numel(idx) < 2, continue; end

        trace = F(:, idx);
        if size(trace,2) < nBins
            trace(:, end+1:nBins) = NaN;
        elseif size(trace,2) > nBins
            trace = trace(:,1:nBins);
        end

        X(:,:,k) = trace;
        usedTrialIdx(k) = true;
    end

    X = X(:,:,usedTrialIdx);
    y = y(usedTrialIdx);

    if ~isempty(X)
        reshaped = reshape(X, nC, []);
        mu = mean(reshaped, 2, 'omitnan');
        X = bsxfun(@minus, X, mu);
    end
end
