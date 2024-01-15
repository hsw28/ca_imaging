function CS_US_id_struct = BULKconverttoframe(US_time_structure, Ca_timestamps)

%converts from a timestamp to a frame #.
%then converts to a spike train (can uncomment this) putting a 1-5 for CS and a 6-10 for US, and a 0 for neither


ratios = struct();
fields_US = fieldnames(US_time_structure);
%fields_CS = fieldnames(CS_time_structure);
fields_CA = fieldnames(Ca_timestamps);

if numel(fields_US) ~= numel(fields_CA)
error('your timestamp and US structures do not have the same number of values. you may need to pad your US structure for exploration days')
end


for i = 1:numel(fields_US)
      field_name = fields_US{i};
      US_timestoconvert = US_time_structure.(field_name);

  %    field_name = fields_CS{i};
  %    CS_timestoconvert = CS_time_structure.(field_name);

      field_name = fields_CA{i};
      Ca_ts = Ca_timestamps.(field_name);

      index = strfind(field_name, '_');
      date = field_name(index(2)+1:end)


            timestamps = Ca_ts;
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

            %  currconv_CS = CS_timestoconvert(k);
            %  [c index] = min(abs(timestamps-currconv_CS));
            %  if (currconv_CS-timestamps(index))>0
            %    CS_frame = ceil(index./2);
            %  else
            %    CS_frame = floor(index./2);
            %  end
                %allframes(CS_frame+1:US_frame-1)=10;
                %allframes(US_frame+0:US_frame+3)=20; %1.7

                %allframes(US_frame-5:US_frame-1)=10;
                %allframes(US_frame+0:US_frame+4)=20;

                if US_frame-5>0 && US_frame+4<=length(allframes)
                  allframes(US_frame-10:US_frame-6)=[-1,-1,-1,-1,-1];
                  allframes(US_frame-5:US_frame-1)=[1,2,3,4,5];
                  allframes(US_frame+0:US_frame+4)=[6,7,8,9,10];
                elseif US_frame-5<=0
                  startpoint = 5+(US_frame-5);
                  allframes(1:startpoint) = [5-startpoint+1:1:5];
                elseif US_frame+4>length(allframes)
                  endpoint = length(allframes)-US_frame + 6;
                  allframes(US_frame+0:end)=[6:1:endpoint];
                end


            end

            tsindex = 2:2:length(timestamps);
            timestamps = timestamps(tsindex);

            if length(timestamps)~=length(allframes)
              warning('your timestamps and frames arent same length')
            end


        CS_US_id_struct.(sprintf('CSUS_id_%s', date)) = [allframes' timestamps]';

  end
