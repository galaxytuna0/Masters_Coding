function [reconstructed_signals, t] = zeroFilt(data, Fs, multisine_frequencies, window_duration_s, plotc)
% zeroFilt Filters a signal by isolating multisine frequencies using FFT/IFFT.
%
%   [reconstructed_signals, t] = zeroFilt(data, Fs, multisine_frequencies)
%   processes a signal 'data' by breaking it into windows.
%   It then uses the FFT and IFFT to isolate and reconstruct individual
%   multisine frequencies specified in 'multisine_frequencies'. This
%   function trims the original data to a length that is a multiple of the
%   window size, discarding any remaining samples.
%
%   Inputs:
%     data                   - A column vector (n x 1) of time-domain data.
%     Fs                     - The sampling frequency in Hz.
%     multisine_frequencies  - A list of frequencies to filter for,
%                              e.g., [0.3, 0.5, 0.7, 1.1].
%
%   Outputs:
%     reconstructed_signals  - A cell array where each cell contains the
%                              reconstructed filtered signal for a
%                              corresponding frequency in 'multisine_frequencies'.
%     t                      - The time vector for the reconstructed signals.

%% Problem 1: Windowing a signal of varying length

% Calculate the total number of samples and the window size
total_samples = numel(data);
window_samples = round(window_duration_s * Fs);

% Calculate the number of full windows and trim the data.
num_full_windows = floor(total_samples / window_samples);
samples_trimmed_data = num_full_windows * window_samples;
data = data(1:samples_trimmed_data);

% Initialize a cell array to store the data from each full window.
data_windows = cell(1, num_full_windows);

% Segment the data into full windows.
fprintf('Processing data in windows...\n');
for i = 1:num_full_windows
    start_index = (i - 1) * window_samples + 1;
    end_index = i * window_samples;
    data_windows{i} = data(start_index:end_index);
    fprintf('  - Extracted Full Window %d (Samples %d to %d)\n', i, start_index, end_index);
end

%% Problem 2: Filtering each frequency using FFT/IFFT with zero-padding

% Define the desired FFT length to maintain a consistent frequency resolution.
desired_fft_length = window_samples;

% Calculate the desired frequency resolution.
desired_freq_res = Fs / desired_fft_length;
fprintf('\nDesired frequency resolution for filtering: %.2f Hz\n', desired_freq_res);

% Initialize a cell array to store the reconstructed signals for each frequency.
reconstructed_signals = cell(1, length(multisine_frequencies));
for j = 1:length(multisine_frequencies)
    reconstructed_signals{j} = []; % Start with an empty array.
end

% Create the time vector for the trimmed data.
t = (0:samples_trimmed_data-1) / Fs;

% Process each window.
for i = 1:num_full_windows
    current_window_data = data_windows{i};
    N_current = length(current_window_data);
    
    % window = hanning(N_current);
    % windowed_data = current_window_data .* window;
    windowed_data = current_window_data;
    
    % Zero-padding is not needed -- processing full windows.
    padded_data = windowed_data;
    
    Y = fft(padded_data);
    f_fft = (0:desired_fft_length-1) * desired_freq_res;
    
    % Filter and append the result to the reconstructed signals.
    for j = 1:length(multisine_frequencies)
        target_freq = multisine_frequencies(j);
        
        [~, pos_idx] = min(abs(f_fft - target_freq));
        [~, neg_idx] = min(abs(f_fft - (Fs - target_freq)));
        
        Y_filtered = zeros(size(Y));
        Y_filtered(pos_idx) = Y(pos_idx);
        Y_filtered(neg_idx) = Y(neg_idx);
        
        y_filtered = ifft(Y_filtered);
        
        % Append the filtered signal to the corresponding reconstructed signal array.
        reconstructed_signals{j} = vertcat(reconstructed_signals{j}, real(y_filtered));
    end
end


%% Visualize the original and reconstructed signals
if plotc == "yes"
    fprintf('\nReconstruction complete. Plotting results...\n');
    figure;
    % Plot the original trimmed signal in the first subplot.
    subplot(length(multisine_frequencies)+1, 1, 1);
    plot(t, data);
    title('Original Trimmed Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    
    % Loop through and plot each reconstructed frequency component.
    for j = 1:length(multisine_frequencies)
        subplot(length(multisine_frequencies)+1, 1, j+1);
        plot(t, reconstructed_signals{j});
        title(sprintf('Reconstructed Component at %.1f Hz', multisine_frequencies(j)));
        xlabel('Time (s)');
        ylabel('Amplitude');
        grid on;
    end
    
    sgtitle('Filtering and Reconstructing the Entire Signal', 'FontSize', 16);

else
    fprintf('\nReconstruction complete. Suppressed plot results...\n');

end