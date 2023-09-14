function allframes = converttoframe(CS_timestoconvert, US_timestoconvert, Ca_timestamps)
%converts from a timestamp to a frame #.
%then converts to a spike train (can uncomment this) putting a 10 for CS and a 20 for US

timestamps = Ca_timestamps;
if isa(timestamps,'table')
  timestamps = table2array(timestamps);
  timestamps = timestamps(:,2);
end

if size(timestamps,2)==3
  timestamps = timestamps(:,2);
end

if timestamps(5)>2
timestamps = timestamps./1000;
end


allframes = zeros(1,floor(length(timestamps)./2));
for k=1:length(US_timestoconvert)

  currconv_US = US_timestoconvert(k);
  [c index] = min(abs(timestamps-currconv_US));
  if (currconv_US-timestamps(index))>0
    US_frame = ceil(index./2);
  else
    US_frame = floor(index./2);
  end

  currconv_CS = CS_timestoconvert(k);
  [c index] = min(abs(timestamps-currconv_CS));
  if (currconv_CS-timestamps(index))>0
    CS_frame = ceil(index./2);
  else
    CS_frame = floor(index./2);
  end
    %allframes(CS_frame:CS_frame+3)=10;
    %allframes(US_frame:US_frame+4)=20;
    allframes(CS_frame:US_frame-1)=10;
    allframes(US_frame+0:US_frame+2)=20;

end
