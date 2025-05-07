function dateList = autoDateList(animal)
    getDate = @(s) regexp(s,'\d{4}_\d{2}_\d{2}','match','once');
    fNames  = fieldnames(animal.Ca_traces);
    dateList = unique(cellfun(getDate, fNames, 'uni',0));
    dateList = sort(dateList);
end
