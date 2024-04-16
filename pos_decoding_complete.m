function [model values error] = pos_decoding_complete(pos_for_model, spikes_for_model, pos_for_decoding, spikes_for_decoding, tdecode, dim, velthreshold, alignmentdata1, alignmentdata2)

both = find(alignmentdata1>0 & alignmentdata2>0);
want1 = (alignmentdata1(both));
want2 = (alignmentdata2(both));
spikes1 = spikes_for_model(want1,:);
spikes2 = spikes_for_decoding(want2,:);







model = pos_decoding_model(pos_for_model, spikes1, tdecode, dim, velthreshold);
[values error] = pos_decoding_with_model(model, spikes2, pos_for_decoding);
