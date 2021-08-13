function f = getspikesignals(unsortedstructure, sortedstructure, timestamps)
%gives you good spike signal trains from extractedSignalsEst and trims length to timestamps

spikepeaks = unsortedstructure.cnmfeAnalysisOutput.extractedSignalsEst;

good = sortedstructure.validCNMFE; %sorted cells
temp = find(good==1);
spikepeaks = spikepeaks(temp, :);


last = timestamps(end)./1000;
last = floor(last*7.5);
if last<length(train_peak_matrix);
  spikepeaks = spikepeaks(:,1:last);
end
f = spikepeaks;
