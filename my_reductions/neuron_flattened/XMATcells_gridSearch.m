function results = XMATcells_gridSearch(animalName, win_vals, Fs_vals, kernels, preprocessingModes, classifiers, balanceVals)
  %% example call
  % results = XMATcells_gridSearch('rat0314', {[0, 0.75], [0.2, 1]}, [3, 5, 7.5], {'linear', 'rbf', 'polynomial'}, {'zscore', 'normalize', 'none'}, {'svm', 'nb', 'tree'}, [true, false]);

if nargin < 1, error('Must provide animalName'); end
if nargin < 2, win_vals = {[0, 0.75]}; end
if nargin < 3, Fs_vals = [7.5]; end
if nargin < 4, kernels = {'linear'}; end
if nargin < 5, preprocessingModes = {'zscore'}; end
if nargin < 6, classifiers = {'svm'}; end
if nargin < 7, balanceVals = [false]; end

animal = evalin('base', animalName);
dateList = autoDateList(animal);
G = animal.alignmentALL;
nDays = numel(dateList);
nAligned = size(G, 2);
if nAligned < nDays
    dateList = dateList(1:nAligned);
    nDays = nAligned;
end

bandMask = abs((1:nDays)' - (1:nDays)) == 1;
lateDays = (nDays-9):nDays;
lateBand = false(nDays,nDays);
lateBand(lateDays, lateDays) = true;
lateBand = lateBand & bandMask;

total = numel(win_vals)*numel(Fs_vals)*numel(kernels)*numel(preprocessingModes)*numel(classifiers)*numel(balanceVals);
count = 0;

results = struct('win', {}, 'Fs', {}, 'kernel', {}, 'prep', {}, 'classifier', {}, 'balanced', {}, ...
    'meanDiag', {}, 'medDiag', {}, 'meanBand', {}, 'medBand', {}, ...
    'meanLateBand', {}, 'medLateBand', {}, 'meanAll', {}, 'medAll', {});

for wi = 1:numel(win_vals)
    win = win_vals{wi};
    for fi = 1:numel(Fs_vals)
        Fs = Fs_vals(fi);
        for ki = 1:numel(kernels)
            kernel = kernels{ki};
            for pi = 1:numel(preprocessingModes)
                prep = preprocessingModes{pi};
                for ci = 1:numel(classifiers)
                    clf = classifiers{ci};
                    for bi = 1:numel(balanceVals)
                        balanceClasses = balanceVals(bi);

                        count = count + 1;
                        fprintf('Running %d of %d: win=[%.2f %.2f], Fs=%.1f, kernel=%s, prep=%s, clf=%s, balanced=%d\n', ...
                            count, total, win(1), win(2), Fs, kernel, prep, clf, balanceClasses);

                        [accMatrixXmat, pMatrixXmat] = ...
                            XMATcells_crossDayDecoder(animalName, win, Fs, kernel, prep, clf, balanceClasses);

                        maskDiag = eye(nDays)==1;

                        results(count).win = win;
                        results(count).Fs = Fs;
                        results(count).kernel = kernel;
                        results(count).prep = prep;
                        results(count).classifier = clf;
                        results(count).balanced = balanceClasses;
                        display(nanmean(pMatrixXmat(maskDiag)))
                        results(count).meanDiag = nanmean(pMatrixXmat(maskDiag));

                        display(nanmedian(pMatrixXmat(maskDiag)))
                        results(count).medDiag = nanmedian(pMatrixXmat(maskDiag));

                        display(nanmean(pMatrixXmat(bandMask)))
                        results(count).meanBand = nanmean(pMatrixXmat(bandMask));

                        display(nanmedian(pMatrixXmat(bandMask)))
                        results(count).medBand = nanmedian(pMatrixXmat(bandMask));

                        display(nanmean(pMatrixXmat(lateBand)))
                        results(count).meanLateBand = nanmean(pMatrixXmat(lateBand));

                        display(nanmedian(pMatrixXmat(lateBand)))
                        results(count).medLateBand = nanmedian(pMatrixXmat(lateBand));

                        display(nanmean(pMatrixXmat(:)))
                        results(count).meanAll = nanmean(pMatrixXmat(:));

                        display(nanmedian(pMatrixXmat(:)))
                        results(count).medAll = nanmedian(pMatrixXmat(:));
                    end
                end
            end
        end
    end
end

results = struct2table(results);
end
