function results = pos_decoding_complete_shuff_trace(pos_for_model, spikes_for_model, CA_timestamps_model, pos_for_decoding, spikes_for_decoding, CA_timestamps_decoding, tdecode, dim, velthreshold, alignmentdata1, alignmentdata2, iters)

both = find(alignmentdata1>0 & alignmentdata2>0);
want1 = (alignmentdata1(both));
want2 = (alignmentdata2(both));
spikes1 = spikes_for_model(want1,:);
spikes2 = spikes_for_decoding(want2,:);

fprintf('fixing positions') %%%
if (pos_for_model(1,1)-pos_for_model(end,1))./length(pos_for_model) < 1
  pos_for_model = convertpostoframe(pos_for_model, CA_timestamps_model);
end


fprintf('fixing positions') %%%
if (pos_for_decoding(1,1)-pos_for_decoding(end,1))./length(pos_for_decoding) < 1
  pos_for_decoding = convertpostoframe(pos_for_decoding, CA_timestamps_decoding);
end



model = pos_decoding_model(pos_for_model, spikes1, tdecode, dim, velthreshold);



[values error] = pos_decoding_with_model(model, spikes2, pos_for_decoding);

results = NaN(2,iters);
for k = 1:iters
    k
    pos_for_decoding = smoothpos(pos_for_decoding);
    vel = ca_velocity(pos_for_decoding);
    want = find(vel(1,:)>=velthreshold);
    pos_want = pos_for_decoding(want,2:3);
    pos_want = shuffle2columns(pos_want);
    pos_for_decoding(want, 2:3) = pos_want;

    dec_error = ca_decodederror(values, pos_for_decoding, tdecode);
    results(1,k) = nanmean(dec_error(:,6));
    results(2,k) = nanmedian(dec_error(:,6));
end

results = results';
