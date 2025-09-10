function out = getAllDayKeys(rat)
    buckets = {'CS_times','US_times','Ca_peaks','Ca_ts','Ca_traces','pos','csus15','csus30','csus45','csus60','csus90'};
    out = {};
    for b = 1:numel(buckets)
        if isfield(rat, buckets{b}) && isstruct(rat.(buckets{b}))
            out = [out; fieldnames(rat.(buckets{b}))]; %#ok<AGROW>
        end
    end
    out = unique(out);
end
