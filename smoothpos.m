function smoothed_positions = smoothpos(positions)

if size(positions,1)<size(positions,2)
  positions = positions';
end

% Extract x and y coordinates
x_positions = positions(:, 2);
y_positions = positions(:, 3);

% Define the standard deviation of the Gaussian filter in cm
sigma = 2; % Adjust as needed

% Apply the Gaussian filter to the x and y coordinates
smoothed_x_positions = imgaussfilt(x_positions, sigma);
smoothed_y_positions = imgaussfilt(y_positions, sigma);

% Combine the smoothed coordinates with the original time column
smoothed_positions = [positions(:, 1), smoothed_x_positions, smoothed_y_positions];
