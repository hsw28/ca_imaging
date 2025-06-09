function results = CrossDayDecoding_UMAP_gridSearch(animalName, win, neighborVals, minDistVals, metricVals)
%trial-flattened 
%example call
% results = CrossDayDecoding_UMAP_gridSearch('rat0314', [0 .75], [5 10 15 20 25 30 40], [0.05 0.1 0.2 0.3 0.4 0.5 0.7], {'euclidean','cosine','manhattan'});

warning('off','all')

    % Preallocate results struct array
    count = 0;
    totalCombos = length(neighborVals) * length(minDistVals) * length(metricVals);

    % Loop over all param combos
    for nIdx = 1:length(neighborVals)
        for mIdx = 1:length(minDistVals)
            for metIdx = 1:length(metricVals)

                count = count + 1;
                fprintf('\n=== Running %d of %d: neighbors=%d, min_dist=%.2f, metric=%s ===\n', ...
                    count, totalCombos, neighborVals(nIdx), minDistVals(mIdx), metricVals{metIdx});

                % TEMPORARILY override run_umap with current params
                % You must edit UMAPtrial_crossDayDecoder to take these params — OR
                % temporarily change run_umap defaults before calling

                % Option 1 (simplest for now):
                % Just edit UMAPtrial_crossDayDecoder to take (neighbors, min_dist, metric)
                % and pass to run_umap

                % For now I assume you will pass these into UMAPtrial_crossDayDecoder:
                [~, accUMAP] = UMAPtrial_crossDayDecoder(animalName, win, false, false, 0, ...
                    neighborVals(nIdx), minDistVals(mIdx), metricVals{metIdx});

                % Now compute metrics:

                nDays = size(accUMAP,1);

                % Off-diagonal mask
                offDiagMask = ~eye(nDays);

                % +/-1 band mask
                bandMask = abs((1:nDays)' - (1:nDays)) <= 1;

                % Late days:
                lateDays = (nDays-2):nDays;
                lateMask = false(nDays,nDays);
                lateMask(lateDays,lateDays) = true;

                % Band in late days:
                lateBandMask = bandMask & lateMask;

                % Metrics:
                meanFull = nanmean(accUMAP(:));
                medFull  = nanmedian(accUMAP(:));

                meanOffDiag = nanmean(accUMAP(offDiagMask));
                medOffDiag  = nanmedian(accUMAP(offDiagMask));

                meanBand = nanmean(accUMAP(bandMask));
                medBand  = nanmedian(accUMAP(bandMask));

                meanLateBand = nanmean(accUMAP(lateBandMask));
                medLateBand  = nanmedian(accUMAP(lateBandMask));

                % Store in struct array:
                s = struct();
                s.n_neighbors = neighborVals(nIdx);
                s.min_dist = minDistVals(mIdx);
                s.metric = metricVals{metIdx};
                s.meanFull = meanFull;
                s.medFull = medFull;
                s.meanOffDiag = meanOffDiag;
                s.medOffDiag = medOffDiag;
                s.meanBand = meanBand;
                s.medBand = medBand;
                s.meanLateBand = meanLateBand;
                s.medLateBand = medLateBand;

                results(count) = s;

            end
        end
    end

    results = struct2table(results);
    fprintf('\n=== Grid search complete! %d combinations tested. ===\n', totalCombos);

end
