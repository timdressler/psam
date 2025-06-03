function tid_psam_plot_flagged_ICs_TD(EEG, sgtitlestr, varargin)
% tid_psam_plot_flagged_ICs_TD - Plots flagged ICs
%
% Usage:
%   tid_psam_plot_flagged_ICs_TD(EEG, sgtitlestr)
%   tid_psam_plot_flagged_ICs_TD(EEG, sgtitlestr, 'SavePath', 'path/to/save.png')
%   tid_psam_plot_flagged_ICs_TD(EEG, sgtitlestr, 'PlotOn', false)
%
% Inputs:
%   EEG - an EEGLAB structure.
%                 EEG must contain information regarding flagged ICs.
%   sgtitlestr - a string used for the title of the plot.
%
% Optional Inputs:
%   'SavePath' - Full file path (including name) to save the plot (optional).
%   'PlotOn' - Logical flag to control whether the plot is shown (default: true).
%
% Outputs:
%   - (None returned)
%
% Plots:
%   Plots all extracted ICs and marks the flagged ones based on EEG.reject.gcompreject.
%
% Tim Dressler, 24.05.2025 (Updated)

% Parse optional inputs
p = inputParser;
addParameter(p, 'SavePath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'PlotOn', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
savePath = p.Results.SavePath;
plotOn = p.Results.PlotOn;

% Check for flagged ICs
if ~isfield(EEG.reject, 'gcompreject') || isempty(EEG.reject.gcompreject)
    error('No information regarding flagged ICs found.');
end

% Get flagged IC indices
flagged_comps = find(EEG.reject.gcompreject);

% Create figure (visible or not depending on PlotOn)
if plotOn
    fig = figure('Units', 'normalized', 'Position', [0 0.0500 1 0.8771]);
else
    fig = figure('Visible', 'off', 'Units', 'normalized', 'Position', [0 0.0500 1 0.8771]);
end

% Plot each IC
numICs = size(EEG.icaweights, 1);
nRows = ceil(sqrt(numICs));
nCols = ceil(numICs / nRows);

for i = 1:numICs
    subplot(nRows, nCols, i);
    topoplot(EEG.icawinv(:, i), EEG.chanlocs, 'numcontour', 0);
    titlestr = sprintf('IC %d', i);
    if ismember(i, flagged_comps)
        titlestr = [titlestr, ' FLAGGED'];
    end
    title(titlestr);
end

sgtitle(sgtitlestr);

% Save figure if a path is specified
if ~isempty(savePath)
    saveas(fig, savePath);
    fprintf('Plot saved to: %s\n', savePath);
end

% Close figure if plot is suppressed
if ~plotOn
    close(fig);
end

end
