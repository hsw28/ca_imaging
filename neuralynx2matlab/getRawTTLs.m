
function events = getRawTTLs(filename)

FieldSelection(1) = 1;%timestamps
FieldSelection(2) = 0;
FieldSelection(3) = 1;%ttls
FieldSelection(4) = 0;
FieldSelection(5) = 0;
ExtractHeader = 1;
ExtractMode = 1;
%ModeArray(1)=fromInd;
%ModeArray(2)=toInd;

fprintf('if you get an error about incorrect outputs your file is prob empty :-p')
[timestamps, ttls, header] = Nlx2MatEV_v3(filename, FieldSelection, ExtractHeader, ExtractMode);
%timestamps = Nlx2MatEV_v3(filename, FieldSelection, ExtractHeader, ExtractMode);

%header


events=zeros(size(ttls,2),2);
events(:,1) = timestamps';
events(:,2) = ttls';
