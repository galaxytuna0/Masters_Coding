% Author: Matt Galassi
% 64 Channel HDsEMG Laplacian Grid Filter - in progress
% Last edited: Dec. 8th 2025
%
% Designed for 8-8L TMSi Grid with REFA6 EEG amp. Re-maps the REFA channels 
% to the TMSi orientation, and runs an average with the adjacent leads
%

%%

d = readmatrix("Session2Test1.csv");

for i = 1:133
    temp_data = d(:,i);
    offset = temp_data - mean(temp_data);
    hp = HighPassFilter(offset,4000,2,15);
    d_filt(:,i) = (hp)';
end

grid1 = [
17 16 15 14 13  9  5  1;
22 21 20 19 18 10  6  2;
27 26 25 24 23 11  7  3;
32 31 30 29 28 12  8  4
];

grid2 = [
33 34 35 36 37 53 57 61;
38 39 40 41 42 54 58 62;
43 44 45 46 47 55 59 63;
48 49 50 51 52 56 60 64
];

% Example data: [64 channels x time]
data = d_filt(:,1:64); % replace with your EMG

% Apply filter
grid1_filt = laplacian_filter_custom(data, grid1);
grid2_filt = laplacian_filter_custom(data, grid2);


%%

function filtered = laplacian_filter_custom(data, grid)
% data: [time_samples x num_channels] 
% grid: [rows x cols] matrix of channel numbers
%
% Returns filtered: [time_samples x num_channels] after subtracting the
% average of adjacent electrodes

[R, C] = size(grid);
[num_samples, num_channels] = size(data);
filtered = zeros(num_samples, num_channels);

% Loop over each grid position
for r = 1:R
    for c = 1:C
        chan = grid(r,c);
        if chan == 0
            continue; % skip empty positions
        end
        
        neighbors = [];
        
        % Check 4 neighbors (up, down, left, right)
        if r > 1 && grid(r-1,c) ~= 0
            neighbors = [neighbors, data(:, grid(r-1,c))];
        end
        if r < R && grid(r+1,c) ~= 0
            neighbors = [neighbors, data(:, grid(r+1,c))];
        end
        if c > 1 && grid(r,c-1) ~= 0
            neighbors = [neighbors, data(:, grid(r,c-1))];
        end
        if c < C && grid(r,c+1) ~= 0
            neighbors = [neighbors, data(:, grid(r,c+1))];
        end
        
        % Subtract the average of neighbors
        if ~isempty(neighbors)
            filtered(:, chan) = data(:, chan) - mean(neighbors, 2);
        else
            filtered(:, chan) = data(:, chan); % no neighbors
        end
    end
end
end
