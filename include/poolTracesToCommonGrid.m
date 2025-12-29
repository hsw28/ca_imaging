function [t_common, emg_mu, emg_sp, vel_mu, vel_sp] = poolTracesToCommonGrid(tMat, emgMat, velMat, spreadType)
% Aligns each rat's mean trace onto a shared time vector and computes mean + spread across rats.
% spreadType: 'std' -> STD across rats
%             'sem' -> SEM across rats

    % pick a reference timebase (assumes all are identical in practice; but we guard anyway)
    % Use the longest, or the first non-empty
    idxRef = find(~cellfun(@isempty,tMat), 1, 'first');
    if isempty(idxRef)
        t_common = []; emg_mu=[]; emg_sp=[]; vel_mu=[]; vel_sp=[];
        return;
    end
    t_common = tMat{idxRef};

    nR = numel(tMat);
    nT = numel(t_common);

    EMG = nan(nR, nT);
    VEL = nan(nR, nT);

    for i = 1:nR
        if isempty(tMat{i}) || isempty(emgMat{i}) || isempty(velMat{i})
            continue;
        end
        ti = tMat{i};

        % If timebases match exactly, direct assign; else interpolate
        if numel(ti)==nT && all(abs(ti - t_common) < 1e-9)
            EMG(i,:) = emgMat{i};
            VEL(i,:) = velMat{i};
        else
            EMG(i,:) = interp1(ti, emgMat{i}, t_common, 'linear', NaN);
            VEL(i,:) = interp1(ti, velMat{i}, t_common, 'linear', NaN);
        end
    end

    % across-rat mean
    emg_mu = nanmean(EMG, 1);
    vel_mu = nanmean(VEL, 1);

    % across-rat spread
    switch lower(spreadType)
        case 'std'
            emg_sp = nanstd(EMG, 0, 1);
            vel_sp = nanstd(VEL, 0, 1);
        case 'sem'
            nEffE = sum(~isnan(EMG),1);
            nEffV = sum(~isnan(VEL),1);
            emg_sp = nanstd(EMG, 0, 1) ./ sqrt(max(1,nEffE));
            vel_sp = nanstd(VEL, 0, 1) ./ sqrt(max(1,nEffV));
        otherwise
            error('spreadType must be "std" or "sem".');
    end
end
