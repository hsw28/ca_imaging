function f = getspikepeaks(unsortedstructure, sortedstructure)
%gives you only sorted spike peaks
%aftrwards you need to still convert spikepeaks to times: spiketimes = converttrain(spikepeaks, length_in_seconds);

%spikepeaks = unsortedstructure.cnmfeAnalysisOutput.extractedPeaks;
spikepeaks = unsortedstructure.cnmfeAnalysisOutput.extractedSignalsEst;

good = sortedstructure.validCNMFE; %sorted cells
temp = find(good==1);
spikepeaks = spikepeaks(temp, :);
%spikepeaks = good;

f = spikepeaks;
