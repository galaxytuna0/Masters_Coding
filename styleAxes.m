% Author: Matt Galassi
% Date: September 28, 2025
% Last Edited: September 28th, 2025
%
% Function formats plot axes to fit with general course expectations.
%
%   Usage: 
%       Call after your full plot code to adjust axes styles
%       Has an option for Legend Title input.
%


function styleAxes(legtitle)
    ax = gca;
    ax.Box = 'off';            % Remove surrounding box
    ax.FontSize = 20;          % Larger font for readability
    ax.LineWidth = 1.5;        % Thicker axis lines
    ax.TickLength = [0.005 0.005]; % Shorter tick marks
    grid off;                  % No grid in 2D plots
    if nargin >= 1
        h = legend('location','northeastoutside', 'Box', 'off');
        title(h, legtitle);
    else
        legend('location','northeastoutside', 'Box', 'off')
    end
end