function tid_psam_plot_rms_bins_TD(EEG, titlestr, varargin)
% tid_psam_plot_rms_bins_TD.m
% 
% Plots RMS of each EEG channel in specified time bins.
%
% Usage:
%   tid_psam_plot_rms_bins_TD(EEG, titlestr)
%   tid_psam_plot_rms_bins_TD(EEG, titlestr, 'BinDuration', 10)
%   tid_psam_plot_rms_bins_TD(EEG, titlestr, 'SavePath', 'path/to/plot.png')
%   tid_psam_plot_rms_bins_TD(EEG, titlestr, 'PlotOn', false)
%
% Inputs:
%   EEG - EEG dataset structure from EEGLAB (containing EEG.data and EEG.srate)
%   titlestr - a string used for the title of the plot
%
% Optional Inputs:
%   'BinDuration' - Bin duration in seconds (default: 10 seconds)
%   'SavePath' - Full file path (including name) to save the plot (optional)
%   'PlotOn' - Logical flag to control whether the plot is shown (default: true)
%
% Description:
%   This function divides continuous EEG data into specified time bins,
%   computes the root mean square (RMS) of each channel within each bin,
%   and optionally displays the result using an imagesc plot. X-axis corresponds to
%   time bins, and Y-axis to EEG channels. Color indicates RMS amplitude.
%   The plot can also be saved to a specified file path.
%
% Tim Dressler, 03.06.25 

% Parse input arguments
p = inputParser;
addParameter(p, 'BinDuration', 10, @(x) isnumeric(x) && x > 0);
addParameter(p, 'SavePath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'PlotOn', true, @(x) islogical(x));
parse(p, varargin{:});

binDurationSec = p.Results.BinDuration;
savePath = p.Results.SavePath;
plotOn = p.Results.PlotOn;

% Check that EEG contains data and sampling rate
if ~isfield(EEG, 'data') || ~isfield(EEG, 'srate')
    error('EEG must contain fields EEG.data and EEG.srate');
end

% Extract data and sampling rate
[data, srate] = deal(EEG.data, EEG.srate);
[numChannels, numSamples] = size(data);

% Define bin size in samples
binSize = binDurationSec * srate;

% Calculate number of full bins
numBins = floor(numSamples / binSize);

% Preallocate RMS matrix (channels x bins)
rmsMatrix = zeros(numChannels, numBins);

% Compute RMS for each bin
for b = 1:numBins
    idxStart = (b-1)*binSize + 1;
    idxEnd = b*binSize;
    segment = data(:, idxStart:idxEnd);
    rmsMatrix(:, b) = sqrt(mean(segment.^2, 2));
end

% Determine channel labels
if isfield(EEG, 'chanlocs') && ~isempty([EEG.chanlocs.labels])
    chanLabels = {EEG.chanlocs.labels};
else
    chanLabels = arrayfun(@num2str, 1:numChannels, 'UniformOutput', false);
    warning('No channel labels found. Using default numeric labels instead.')
end

% Plot the RMS matrix if requested
if plotOn
    figure;
    imagesc(rmsMatrix);
    colorbar;
    xlabel(sprintf('%.0f s Bins', binDurationSec));
    ylabel('Channels');
    title(titlestr);
    set(gca, 'YDir', 'normal');
    set(gca, 'YTick', 1:numChannels, 'YTickLabel', chanLabels);

    % Save plot if path is provided
    if ~isempty(savePath)
        saveas(gcf, savePath);
        fprintf('Plot saved to: %s\n', savePath);
    end
else
    % Save plot without displaying if needed
    if ~isempty(savePath)
        f = figure('Visible', 'off');
        imagesc(rmsMatrix);
        colorbar;
        xlabel(sprintf('%.0f s Bins', binDurationSec));
        ylabel('Channels');
        title(titlestr);
        set(gca, 'YDir', 'normal');
        set(gca, 'YTick', 1:numChannels, 'YTickLabel', chanLabels);
        saveas(f, savePath);
        close(f);
        fprintf('Plot saved to: %s\n', savePath);
    else 
        disp('Since no savePath is provided and the plot is suppressed nothing is saved and/or shown.')
    end
end

end
