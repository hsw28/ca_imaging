function f = ca_placeevent(events, pos)

% Extract timestamps from position data
postime = pos(:,1);

% Initialize the output matrix for event positions
event_pos = NaN(3,length(events));

% Use knnsearch to find indices of closest timestamps in pos for each event
closestIndices = knnsearch(postime, events');

% Assign the closest position coordinates to each event
event_pos(2:3,:) = pos(closestIndices, 2:3)';  % Transpose to match dimensions
event_pos(1,:) = events;  % Store the original events times

% Return the formatted event positions
f = event_pos;
