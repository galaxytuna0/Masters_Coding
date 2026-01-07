% Author: Matt Galassi
% Date: July 25th, 2025
% Last Edited: September 26th, 2025
%
% (Modified function "convertPoly5ToCsvMultiple" to be executable)
%
% This function allows the user to select multiple .poly5 EEG files
% and converts the EEG data from each file into a separate CSV file.
%
% Need Biosig toolbox for Poly5, or adjust pop_biosig() code for other Poly5 reader

% function convertPoly5ToCsvMultiple()

% Start EEGLAB, but suppress the GUI for scripting purposes
% We start EEGLAB once for the session.
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% Clear any existing datasets in EEGLAB's memory (optional, good practice)
ALLEEG = [];
EEG = [];
CURRENTSET = 0;

% --- 2. Prompt User for .poly5 File Paths ---
% Use uigetfile to allow the user to select multiple .poly5 files
% 'MultiSelect', 'on' enables selection of multiple files.
[filenames, filepath] = uigetfile('*.poly5', 'Select your .poly5 EEG files', 'MultiSelect', 'on');

% Check if the user cancelled the file selection
if isequal(filenames,0) || isequal(filepath,0)
   disp('File selection cancelled. Script aborted.');
   % Clean up EEGLAB environment before exiting
   close all;
   clear ALLEEG EEG CURRENTSET ALLCOM;
   return; % Exit the script
end

% Ensure filenames is a cell array, even if only one file was selected
if ischar(filenames)
    filenames = {filenames};
end

fprintf('Selected %d .poly5 files.\n', length(filenames));

% --- Loop through each selected file ---
for i = 1:length(filenames)
    currentFilename = filenames{i};
    fullPoly5Path = fullfile(filepath, currentFilename);
    fprintf('\n--- Processing file: %s ---\n', fullPoly5Path);

    % Clear EEG structure for the new file to avoid residual data
    EEG = [];

    % --- 3. Load .poly5 Data using pop_biosig ---
    disp('Loading EEG data using pop_biosig...');
    try
        EEG = pop_biosig(fullPoly5Path, 'importevent','off');
        fprintf('Successfully loaded %d channels and %d data points.\n', EEG.nbchan, EEG.pnts);
    catch ME
        warning('Failed to load %s using pop_biosig. Error: %s', currentFilename, ME.message);
        disp('Please ensure the file is not corrupted and EEGLAB biosig plugin is up-to-date. Skipping this file.');
        continue; % Skip to the next file
    end

    % --- 4. Extract EEG Data ---
    eegData = EEG.data;
    eegDataForCsv = eegData'; % Transpose for (time points x channels)
    fprintf('Extracted EEG data dimensions (rows=time points, columns=channels): %s\n', mat2str(size(eegDataForCsv)));

    % --- 5. Determine CSV Output File Name ---
    name = 'quietsitting';
    [~, nameOnly, ~] = fileparts(currentFilename);
    defaultCsvName = [name, '_eeg_data.csv'];
    
    % Use uiputfile for each file to allow user to specify output path/name
    % for individual files, or comment out for automatic saving in source folder
    [csvFilename, csvFilepath] = uiputfile('*.csv', sprintf('Save EEG data for %s as CSV', currentFilename), fullfile(filepath, defaultCsvName));
    
    if isequal(csvFilename,0) || isequal(csvFilepath,0)
       disp('CSV save cancelled for this file. Skipping.');
       continue; % Skip to the next file
    end

    fullCsvPath = fullfile(csvFilepath, csvFilename);
    fprintf('Saving EEG data to: %s\n', fullCsvPath);

    % --- 6. Write EEG Data to CSV File ---
    disp('Writing data to CSV...');
    try
        writematrix(eegDataForCsv, fullCsvPath);
        disp('EEG data successfully written to CSV!');
    catch ME
        warning('Failed to write data to CSV file %s. Error: %s', fullCsvPath, ME.message);
        disp('Please check file permissions and disk space. Skipping this file.');
        continue; % Skip to the next file
    end

end % End of for loop

disp('\nAll selected files processed.');

% --- Clean up EEGLAB environment after all files are processed ---
eeglab redraw; % Refreshes EEGLAB, ensures it's in a clean state
close all;     % Closes all figures, including any EEGLAB main window
clear ALLEEG EEG CURRENTSET ALLCOM; % Clears EEGLAB variables from workspace

disp('Script finished.');

% end