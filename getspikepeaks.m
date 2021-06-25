function f = getspikepeaks(unsortedstructure, sortedstructure)
%gives you only sorted spike peaks

spikepeaks = unsortedstructure.cnmfeAnalysisOutput.extractedPeaks;

good = sortedstructure.validCNMFE; %sorted cells
temp = find(good==1);
spikepeaks = spikepeaks(temp, :);

f = spikepeaks;
