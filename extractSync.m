% Author: Matt Galassi
% Date: September 25, 2025
% Last Edited: October 3rd, 2025
%
% Function extracts the EMG and force plate data from a LabChart data file.
% Mat-file should contain only 1 trial.
%
%   Usage: 
%       xTime, Channel1, Channel2, FP, FCR] = extractSync('filename')
%       Specify the exported .mat file as input.
%
%   Output:
%       xTime: Time base for trial
%       Channel1: placeholder
%       Channel1: placeholder
%       FP: force plate data z-axis
%       FCR: EMG signal for Flexor Carpi Radialis



function [xTime, Channel1, Channel2, FP, FCR] = extractSync(filename)

% Load data file
datafile = load(filename);

% Extract data channels
for ch = 1:4
    chStart(ch) = datafile.datastart(ch);
    chEnd(ch) = datafile.dataend(ch);
end

datasamples = (chEnd(1) - chStart(1) + 1);
xTime = (-0.1*datafile.samplerate(1)):((datasamples-1) - 0.1*datafile.samplerate(1));
xTime = xTime/datafile.samplerate(1);

% Output
Channel1 = datafile.data(chStart(1):chEnd(1));
Channel2 = datafile.data(chStart(2):chEnd(2));
FP = datafile.data(chStart(3):chEnd(3));
FCR = datafile.data(chStart(4):chEnd(4));


end


