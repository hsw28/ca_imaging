function [shuffMeans, shuffSEMs] = plotRateMapSimilaritySummary(nShuffle)
% ------------------------------------------------------------------------
% COMBINED VERSION:
% 1. Logic: Uses your working 'ismembertol' logic for Task/NonTask split.
% 2. Velocity: Filters Non-Task spikes for v >= 4 cm/s.
% 3. Maps: Uses FULL 'pos' for both to ensure dimensions match (fixes NaNs).
% ------------------------------------------------------------------------
if nargin < 1, nShuffle = 500; end


ratNames = {'rat0222', 'rat0307', 'rat0313', 'rat0314', 'rat0816'};
metrics  = {'Pearson', 'Spearman', 'Cosine'};
metrics  = {'Pearson'};



nRats    = numel(ratNames);
nMetrics = numel(metrics);


actualMeans = nan(nRats, nMetrics);
shuffMeans  = nan(nRats, nMetrics);
actualSEMs  = nan(nRats, nMetrics);
shuffSEMs   = nan(nRats, nMetrics);
pRat        = nan(nRats, nMetrics);
fracSig     = nan(nRats, nMetrics);

for r = 1:nRats
  r = r
    try
        rat = evalin('base', ratNames{r});
    catch
        fprintf('Rat %s not found. Skipping.\n', ratNames{r});
        continue;
    end

    dates = autoDateList(rat);
    idx   = find(strcmp(dates, rat.An));
    these = dates(idx-2:idx);

    actualCorrs = [];
    shuffCorrs  = [];
    pNeuron     = [];

    % Shuffle accumulators
    sumRep   = zeros(nShuffle, nMetrics);
    countRep = zeros(nShuffle, nMetrics);

    for d = 1:numel(these)
      d = d
        dateStr = these{d};
        % ---- Load Data ----
        try
            pos = rat.pos.(['pos_' dateStr]);
            cs  = rat.CS_times.(['CS_' dateStr]);
            spkMat = rat.Ca_peaks.(['CA_peaks_' dateStr]);
            mask = rat.ratemask.(['ratemask_' dateStr]);
        catch
            continue;
        end

        % ---- Clean Position ----
        ts = pos(:,1);
        xy = pos(:,2:3);
        % Sort time if needed
        if any(diff(ts) <= 0)
            [ts, ord] = sort(ts);
            xy = xy(ord,:);
            pos = [ts, xy];
        end

        % ---- 1. Define Task Times ( [CS, CS+2] ) ----
        taskMask = false(size(ts));
        cs = cs(isfinite(cs));
        for i = 1:numel(cs)
            taskMask = taskMask | (ts >= cs(i) & ts <= cs(i) + 2);
        end
        taskTS = ts(taskMask); % Timestamps inside task

        % ---- 2. Prepare Velocity (for filtering non-task) ----
        vel = ca_velocity(pos); % [t; v]
        hasVel = ~isempty(vel);
        if hasVel
            vt = vel(1,:)'; vmag = vel(2,:)';
            % Remove duplicates for interp1
            [vt, uIdx] = unique(vt);
            vmag = vmag(uIdx);
        end

        % ---- Neuron Loop ----
        keptNeur = find(mask(:)==1);

        for ni = keptNeur(:)'
            spikes = spkMat(ni,:);
            spikes = spikes(isfinite(spikes));
            if numel(spikes) < 5, continue; end

            % A) Task Spikes: Use your working ismembertol logic
            % (Finds spikes that occur during task timestamps)
            taskSpikes = spikes(ismembertol(spikes, taskTS, 1e-3));

            % B) Non-Task Spikes: Everything else
            nonTaskSpikes_All = setdiff(spikes, taskSpikes, 'stable');

            % C) Velocity Gate the Non-Task Spikes
            if hasVel && ~isempty(nonTaskSpikes_All)
                % Interpolate velocity at spike times
                spkV = interp1(vt, vmag, nonTaskSpikes_All, 'linear', 0);
                nonTaskSpikes = nonTaskSpikes_All(spkV >= 4);
            else
                nonTaskSpikes = nonTaskSpikes_All;
            end

            % Check spike counts
            if numel(taskSpikes) < 3 || numel(nonTaskSpikes) < 3, continue; end

            % ---- Actual Maps (Use FULL pos for both to match dims) ----
            R1 = CA_normalizePosData(taskSpikes(:), pos, 2.5, 1);
            R2 = CA_normalizePosData(nonTaskSpikes(:), pos, 2.5, 1);

            R1 = imgaussfilt(R1, 0.75); R2 = imgaussfilt(R2, 0.75);

            % Flatten
            v1 = R1(:); v2 = R2(:);

            % Remove NaNs and Zeros (Occupancy gaps)
            % Crucial: We only correlate pixels where the rat VISITED in both conditions
            % (or at least visited in the session, depending on normalizePosData internals)
            good = isfinite(v1) & isfinite(v2) & (v1~=0 | v2~=0);

            if nnz(good) < 10, continue; end
            v1 = v1(good); v2 = v2(good);

            % Zero Variance Check (Prevents NaNs in corr)
            if std(v1)==0 || std(v2)==0, continue; end

              act = [corr(v1,v2,'type','Pearson')];
      %      act = [corr(v1,v2,'type','Pearson'), ...
      %             corr(v1,v2,'type','Spearman'), ...
      %             (dot(v1,v2)/(norm(v1)*norm(v2)))];

            if any(isnan(act)), continue; end
            actualCorrs(end+1,:) = act; %#ok<AGROW>


            % ---- Shuffles ----
            allSpk = [taskSpikes(:); nonTaskSpikes(:)];
            lbl    = [ones(numel(taskSpikes),1); zeros(numel(nonTaskSpikes),1)];

            shThis = nan(nShuffle, nMetrics);

            parfor s = 1:nShuffle
                perm = lbl(randperm(numel(lbl)));
                s1 = allSpk(perm==1);
                s2 = allSpk(perm==0);

                % Use FULL pos for shuffles too
                Rm1 = CA_normalizePosData(s1, pos, 2.5, 1);
                Rm2 = CA_normalizePosData(s2, pos, 2.5, 1);

                Rm1 = imgaussfilt(Rm1, 0.75); Rm2 = imgaussfilt(Rm2, 0.75);

                vs1 = Rm1(:); vs2 = Rm2(:);
                gS = isfinite(vs1) & isfinite(vs2) & (vs1~=0 | vs2~=0);

                if nnz(gS) >= 10
                    vs1=vs1(gS); vs2=vs2(gS);
                    if std(vs1)==0 || std(vs2)==0
                        shThis(s,:) = 0;
                    else
                      %  shThis(s,:) = [corr(vs1,vs2,'type','Pearson'), ...
                      %                 corr(vs1,vs2,'type','Spearman'), ...
                      %                 (dot(vs1,vs2)/(norm(vs1)*norm(vs2)))];
                        shThis(s,:) = [corr(vs1,vs2,'type','Pearson')];
                    end
                else
                    shThis(s,:) = 0;
                end
            end

            shMean = nanmean(shThis,1);
            shMean(isnan(shMean)) = 0;
            shuffCorrs(end+1,:) = shMean; %#ok<AGROW>
            pNeuron(end+1,:) = 1 - mean(shThis >= act, 1, 'omitnan'); %#ok<AGROW>

            valid = ~isnan(shMean);
            sumRep(valid) = sumRep(valid) + shMean(valid);
            countRep(valid) = countRep(valid) + 1;
        end
    end

    if isempty(actualCorrs)
        fprintf('No valid neurons for %s\n', ratNames{r});
        continue;
    end

    actualMeans(r,:) = nanmean(actualCorrs,1);
    shuffMeans(r,:)  = nanmean(shuffCorrs,1);
    actualSEMs(r,:)  = nanstd(actualCorrs,0,1)/sqrt(size(actualCorrs,1));
    shuffSEMs(r,:)   = nanstd(shuffCorrs,0,1)/sqrt(size(shuffCorrs,1));

    % Rat-level P-Value
    ratShuffMean = sumRep ./ countRep;
    ratShuffMean(~isfinite(ratShuffMean)) = NaN;


    for m = 1:nMetrics
      pRat(r,m) = length(find(ratShuffMean(:,m)<actualMeans(r,m)))./(length(isfinite(ratShuffMean(:,m))));
    end
    fracSig(r,:) = mean(pNeuron < 0.05, 1, 'omitnan');
end

% ---- Plotting ----
figure('Color','w','Position',[200 300 1200 400]);
for m = 1:nMetrics
    subplot(1,3,m); cla; hold on;

    bData = [actualMeans(:,m), shuffMeans(:,m)];
    errData = [actualSEMs(:,m), shuffSEMs(:,m)];

    if all(isnan(bData(:))), continue; end

    b = bar(bData, 'grouped');
    if numel(b) > 1
        b(1).FaceColor = [0.3 0.5 0.9];
        b(2).FaceColor = [0.7 0.7 0.7];
    end

    % Safe Error Bars
    for k = 1:numel(b)
        xEnd = b(k).XEndPoints;
        errorbar(xEnd, bData(:,k), errData(:,k), 'k.', 'CapSize',0);
    end

    xticks(1:nRats); xticklabels(ratNames);
    ylabel('Similarity'); title(metrics{m}); ylim([0 1.1]);

    pRat
    % Rat-level stars
    for rr = 1:nRats
        val = pRat(rr,m);
        if ~isnan(val)
            if val<0.001, s='***'; elseif val<0.01, s='**'; elseif val<0.05, s='*'; else, s=''; end
            if ~isempty(s)
                h = max(bData(rr,:));
                text(rr, h*1.05, s, 'HorizontalAlignment','center','FontSize',14);
            end
        end
    end
end

% Plot Fraction Significant Neurons
figure('Color','w','Position',[250 320 1200 350],'Name','Fraction sig. neurons');
for m = 1:nMetrics
    subplot(1,3,m);
    bar(fracSig(:,m)*100,'FaceColor',[0.4 0.7 0.4]);
    ylim([0 100]); ylabel('% neurons p<0.05');
    xticks(1:nRats); xticklabels(ratNames);
    title([metrics{m} ' : sig-neuron fraction']);
    yline(5,'--k','\alpha=0.05','LabelHorizontalAlignment','left');
end

fprintf('\nDone.\n');
end
