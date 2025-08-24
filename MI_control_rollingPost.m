function rolling_struct = MI_control_rollingPost(spike_structure, pos_structure, velthreshold, dim, CA_timestamps, CSUS_id_struct, nIter, winEdges)
% Rolling post-CS/US exclusion:
% - Deterministic: remove spikes in the CURRENT post window from the non-CS/US pool
% - Control: remove the same # of spikes from POST-ONLY spikes OUTSIDE the current window
% Labels: 1–7 = CS/US block, 0 = baseline, >7 = post. Seconds = (label-7)/7.5.

if nargin < 7 || isempty(nIter),   nIter   = 100; end
if nargin < 8 || isempty(winEdges), winEdges = [0 2; 1 3; 2 4; 3 5; 4 6; 5 7; 8 10; 9 11; 10 12; 11 13]; end

set(0,'DefaultFigureVisible','off');

fields_spikes = fieldnames(spike_structure);
fields_pos    = fieldnames(pos_structure);
fields_cats   = fieldnames(CA_timestamps);
fields_CSUS   = fieldnames(CSUS_id_struct);

% last 3 days
fields_spikes = fields_spikes(max(1,end-2):end);
fields_pos    = fields_pos(max(1,end-2):end);
fields_cats   = fields_cats(max(1,end-2):end);
fields_CSUS   = fields_CSUS(max(1,end-2):end);

rolling_struct = struct();

for iDay = 1:numel(fields_spikes)
    fieldName_sp   = fields_spikes{iDay}
    fieldName_pos  = fields_pos{iDay};
    fieldName_ts   = fields_cats{iDay};
    fieldName_CSUS = fields_CSUS{iDay};

    peaks_time = spike_structure.(fieldName_sp);   % Ncells x NspikeTimes
    pos        = pos_structure.(fieldName_pos);    % [t x y]
    CSUS_id    = CSUS_id_struct.(fieldName_CSUS); % [labels; times]

    % day tag
    us_ = strfind(fieldName_sp, '_');
    spikes_date = fieldName_sp(us_(2)+1:end);

    % ---- velocity & interpolation ----
    % (dedupe pos first for clean velocity + interpolation)
    [~, uniqPos] = unique(pos(:,1), 'stable');
    pos = pos(uniqPos,:);

    vel      = ca_velocity(pos);         % [speed; time]
    vel_time = vel(2,:)';
    vel_mag  = vel(1,:)';

    % labels onto vel_time
    interp_CSUS = interp1(CSUS_id(2,:), CSUS_id(1,:), vel_time, 'nearest', 0);

    % pos to vel_time
    ix = interp1(pos(:,1), pos(:,2), vel_time, 'linear', NaN);
    iy = interp1(pos(:,1), pos(:,3), vel_time, 'linear', NaN);

    % occupancy mask: movement & non-CS/US (0 or >7)
    is_nonCSUS_occ = (interp_CSUS <= 0 | interp_CSUS > 7.0);
    validIdx = (vel_mag >= velthreshold) & is_nonCSUS_occ & ~isnan(ix) & ~isnan(iy);
    posDat   = [vel_time(validIdx), ix(validIdx), iy(validIdx)];

    nCells = size(peaks_time,1);
    nWins  = size(winEdges,1);

    MI_base = nan(1, nCells);
    MI_win  = nan(nWins, nCells);
    MI_rand = nan(nWins, nCells);

    % ---------- MI_base over the *non-CS/US movement* pool ----------
    for k = 1:nCells
        spk = peaks_time(k,:); spk = spk(~isnan(spk) & spk > 0);
        if isempty(spk), continue; end

        spike_vel = interp1(vel_time, vel_mag, spk, 'linear');
        csus_spk  = interp1(CSUS_id(2,:), CSUS_id(1,:), spk, 'previous', 0);

        nonCSUS_spikes = spk((spike_vel >= velthreshold) & (csus_spk <= 0 | csus_spk > 7.0));
        if isempty(nonCSUS_spikes) || size(posDat,1) < 10, continue; end

        [~,~,~,~,spikeprob,occprob] = CA_normalizePosData(nonCSUS_spikes, posDat, dim, 1.0);
        if size(spikeprob,1) < size(spikeprob,2), spikeprob = spikeprob'; end
        if size(occprob,1)  < size(occprob,2),  occprob  = occprob';  end
        MI_base(k) = mutualinfo([spikeprob, occprob]);
    end

    % ---------- Rolling windows ----------
    parfor k = 1:nCells
        spk = peaks_time(k,:); spk = spk(~isnan(spk) & spk > 0);

        MI_win_k  = nan(nWins,1);
        MI_rand_k = nan(nWins,1);

        if isempty(spk) || isnan(MI_base(k)) || isempty(posDat) || size(posDat,1) < 10
            MI_win(:,k)  = MI_win_k;
            MI_rand(:,k) = MI_rand_k;
            continue;
        end

        spike_vel = interp1(vel_time, vel_mag, spk, 'linear');
        csus_spk  = interp1(CSUS_id(2,:), CSUS_id(1,:), spk, 'previous', 0);

        nonCSUS_spikes = spk((spike_vel >= velthreshold) & (csus_spk <= 0 | csus_spk > 7.0));

        % post mask & seconds since label==7 (ticks/7.5)
        post_mask = (spike_vel >= velthreshold) & (csus_spk > 7.0);
        post_sec  = (csus_spk - 7.0) ./ 7.5;

        for w = 1:nWins
            % ---- define a half-open window [start, end) to avoid double counting ----
            wstart = winEdges(w,1);
            wend   = winEdges(w,2);

            inWin_post = post_mask & (post_sec > wstart) & (post_sec <= wend);
            spk_inWin_all = spk(inWin_post);

            % deterministic pool = non-CS/US movement spikes
            spk_inWin = intersect(spk_inWin_all, nonCSUS_spikes);
            nRemove   = numel(spk_inWin);

            if nRemove <2
                % No info for this cell×window → keep as NaN (prevents biased early bump)
                MI_win_k(w)  = NaN;
                MI_rand_k(w) = NaN;
                continue;
            end

            % ---------- deterministic removal (current window only) ----------
            keep_win = setdiff(nonCSUS_spikes, spk_inWin, 'stable');
            if numel(keep_win) >= 2
                [~,~,~,~,spikeprobW,occprobW] = CA_normalizePosData(keep_win, posDat, dim, 1.0);
                if size(spikeprobW,1) < size(spikeprobW,2), spikeprobW = spikeprobW'; end
                if size(occprobW,1)  < size(occprobW,2),  occprobW  = occprobW';  end
                MI_win_k(w) = mutualinfo([spikeprobW, occprobW]);
            else
                MI_win_k(w) = NaN;
            end

            % ---------- POST-ONLY control, excluding the current window ----------
            post_outWin = post_mask & ~(post_sec > wstart & post_sec <= wend);
            % IMPORTANT: control pool must ALSO live in the non-CS/US movement pool
            post_pool   = intersect(spk(post_outWin), nonCSUS_spikes);

            mi_accum = 0; n_ok = 0;
            if numel(post_pool) >= nRemove
                for it = 1:nIter
                    rand_rm = randsample(post_pool, nRemove, false);
                    keep_rd = setdiff(nonCSUS_spikes, rand_rm, 'stable');
                    if numel(keep_rd) < 2, continue; end
                    [~,~,~,~,spikeprobR,occprobR] = CA_normalizePosData(keep_rd, posDat, dim, 1.0);
                    if size(spikeprobR,1) < size(spikeprobR,2), spikeprobR = spikeprobR'; end
                    if size(occprobR,1)  < size(occprobR,2),  occprobR  = occprobR';  end
                    mi_accum = mi_accum + mutualinfo([spikeprobR, occprobR]);
                    n_ok = n_ok + 1;
                end
            end
            if n_ok > 0
                MI_rand_k(w) = mi_accum / n_ok;
            else
                MI_rand_k(w) = NaN;   % no valid control → leave missing
            end
        end

        MI_win(:,k)  = MI_win_k;
        MI_rand(:,k) = MI_rand_k;
    end


    out.winEdges = winEdges;
    out.MI_win   = MI_win;
    out.MI_rand  = MI_rand;
    out.MI_base  = MI_base;
    out.dMI      = MI_win - MI_rand;


    rolling_struct.(sprintf('MI_%s', spikes_date)) = out;
end
end
