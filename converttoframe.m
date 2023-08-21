function allframes = converttoframe(US_timestoconvert, timestamps)
%converts from a timestamp to a frame #.
%then converts to a spike train (can uncomment this) putting a 10 for CS and a 20 for US

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
  if (currconv_US-c)>0
    US_frame = ceil(index./2);
  else
    US_frame = floor(index./2);
  end

  currconv_CS = currconv_US-.75;
  [c index] = min(abs(timestamps-currconv_CS));
  if (currconv_CS-c)>0
    CS_frame = ceil(index./2);
  else
    CS_frame = floor(index./2);
  end

    allframes(CS_frame:US_frame)=10; %makes a spike train
    allframes(US_frame:US_frame+4)=20;
end
