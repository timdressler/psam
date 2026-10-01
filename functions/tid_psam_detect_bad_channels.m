function [chans_to_interp, plot_handle] = tid_psam_detect_bad_channels(EEG, bemobil_params, varargin)
% tid_psam_detect_bad_channels.m
%
% Detects bad EEG channels using BeMoBIL (clean_rawdata).
% Generates a standardized quality control plot highlighting rejected channels.
%
% Usage:
%   [chans_to_interp, plot_handle] = tid_psam_detect_bad_channels(EEG, [], 'PlotOn', true)
%   [chans_to_interp, plot_handle] = tid_psam_detect_bad_channels(EEG, bemobil_params, 'PlotOn', false)
%   [chans_to_interp, plot_handle] = tid_psam_detect_bad_channels(EEG, [], 'SegmentDuration', 'fulldata')
%   [chans_to_interp, plot_handle] = tid_psam_detect_bad_channels(EEG, [], 'NumSegments', 3, 'SegmentDuration', 10)
%
% Input:
%   EEG            - The EEGLAB EEG structure.
%   bemobil_params - Struct containing parameters for BeMoBIL bad channel
%                    detection. If empty [], defaults are used.
%
% Optional Name-Value Inputs:
%   'PlotOn'          - Boolean. If true, generates and returns a handle to a 
%                       sanity check plot. Default = true.
%   'SegmentDuration' - Duration of each plotted segment in seconds, or 'fulldata' 
%                       to plot the entire dataset. Default = 5.
%   'NumSegments'     - Number of random time segments to plot. Default = 4.
%
% Output:
%   chans_to_interp - Array of channel indices marked for rejection.
%   plot_handle     - Figure handle for the rejection plot (empty if PlotOn=false).
%
% Note. Requires to have bemobil_detect_bad_channels_modified_TD() available!
%
% Tim Dressler, 2026

%% 1. Parse Inputs and Setup Defaults
if nargin < 2
    error('Not enough input arguments. Usage: tid_psam_detect_bad_channels(EEG, bemobil_params, ...)');
end

p = inputParser;
addParameter(p, 'PlotOn', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SegmentDuration', 5, @(x) (isnumeric(x) && x > 0) || ((ischar(x) || isstring(x)) && strcmpi(x, 'fulldata')));
addParameter(p, 'NumSegments', 4, @(x) isnumeric(x) && x > 0);
parse(p, varargin{:});

plotOn = p.Results.PlotOn;
seg_duration = p.Results.SegmentDuration;
n_segments = p.Results.NumSegments;

% Set BeMoBIL Defaults
if isempty(bemobil_params)
    bemobil_params.chancorr_crit = 0.8;
    bemobil_params.chan_max_broken_time = 0.5;
    bemobil_params.chan_detect_num_iter = 10;
    bemobil_params.chan_detected_fraction_threshold = 0.5;
    bemobil_params.num_chan_rej_max_target = 1/5;
    bemobil_params.flatline_crit = 'off';
    bemobil_params.line_noise_crit = 'off';
end

plot_handle = [];

%% 2. Run BeMoBIL Detection
fprintf('Running BeMoBIL Bad Channel Detection...\n');
[bemobil_bad_chans, ~, ~, ~, ~] = bemobil_detect_bad_channels_modified_TD(EEG, [], 1, ...
    bemobil_params.chancorr_crit, bemobil_params.chan_max_broken_time, ...
    bemobil_params.chan_detect_num_iter, bemobil_params.chan_detected_fraction_threshold, ...
    bemobil_params.num_chan_rej_max_target, bemobil_params.flatline_crit, ...
    bemobil_params.line_noise_crit);

%% 3. Sanitize Output
chans_to_interp = double(bemobil_bad_chans(:)');
all_chans = 1:EEG.nbchan;
good_chans = setdiff(all_chans, chans_to_interp);

fprintf('BeMoBIL flagged %d channels for removal.\n', length(chans_to_interp));

%% 4. Detailed All-Channel Plotting
if plotOn
    plot_handle = figure('Visible', 'off', 'Color', 'w', 'Renderer', 'painters');
    set(plot_handle, 'Units', 'pixels', 'Position', [0, 0, 1920, 1080]);
    
    % Define Colors
    col_good = [0.75 0.75 0.75]; % Light Grey
    col_bad = [0.85 0.125 0.048]; % Red
    
    % Setup Title
    sgtitle(plot_handle, sprintf('Bad Channel Detection (BeMoBIL)\nChannels Removed: %d', length(chans_to_interp)), ...
        'Interpreter', 'none', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Time Segment Setup
    srate = EEG.srate;
    is_fulldata = (ischar(seg_duration) || isstring(seg_duration)) && strcmpi(seg_duration, 'fulldata');
    
    if is_fulldata
        n_samples = EEG.pnts;
        n_subplots = 1;
        start_points = 1;
    else
        n_samples = round(seg_duration * srate);
        n_subplots = n_segments;
        
        if EEG.pnts > (n_samples * n_subplots)
            valid_starts = 1:(EEG.pnts - n_samples);
            start_points = sort(valid_starts(randperm(length(valid_starts), n_subplots)));
        else
            start_points = 1;
            n_samples = EEG.pnts;
            n_subplots = 1;
            is_fulldata = true; % Fallback if data is too short
        end
    end
    
    % Channel Labels Fallback
    if isfield(EEG, 'chanlocs') && ~isempty([EEG.chanlocs.labels])
        chanLabels = {EEG.chanlocs.labels}';
    else
        chanLabels = arrayfun(@(x) sprintf('Ch%d', x), 1:EEG.nbchan, 'UniformOutput', false)';
    end
    
    % Global Vertical Offset Calculation
    signal_swing = prctile(EEG.data(:), 95) - prctile(EEG.data(:), 5);
    offset_val = max(signal_swing * 1.2, 30);
    offsets = (EEG.nbchan:-1:1)' * offset_val;
    y_font_size = max(4, min(9, 600 / EEG.nbchan));
    
    for i = 1:n_subplots
        subplot(1, n_subplots, i);
        hold on;
        sp = start_points(i);
        ep = sp + n_samples - 1;
        time_vec = (sp:ep) / srate;
        
        % Plot Dummy Lines for Legend
        if i == 1
            h1 = plot(NaN, NaN, 'Color', col_good, 'LineWidth', 1.5);
            h2 = plot(NaN, NaN, 'Color', col_bad, 'LineWidth', 2.0);
            legend([h1, h2], {'Good (Kept)', 'Bad (Dropped)'}, ...
                'Location', 'northoutside', 'Orientation', 'horizontal');
        end
        
        % Plot each channel
        for c = 1:EEG.nbchan
            trace = EEG.data(c, sp:ep) - mean(EEG.data(c, sp:ep)) + offsets(c);
            if ismember(c, chans_to_interp)
                plot(time_vec, trace, 'Color', col_bad, 'LineWidth', 0.4, 'HandleVisibility', 'off');
            else
                plot(time_vec, trace, 'Color', col_good, 'LineWidth', 0.4, 'HandleVisibility', 'off');
            end
        end
        
        % Formatting
        xlim([time_vec(1), time_vec(end)]);
        ylim([0, (EEG.nbchan + 1) * offset_val]);
        
        if is_fulldata
            title(sprintf('Full Data Timecourse: %ds - %ds', round(time_vec(1)), round(time_vec(end))));
        else
            title(sprintf('Segment: %ds - %ds', round(time_vec(1)), round(time_vec(end))));
        end
        
        xlabel('Time (s)');
        box off;
        
        if i == 1
            yticks(flip(offsets));
            yticklabels(flip(chanLabels));
            set(gca, 'TickLabelInterpreter', 'none', 'FontSize', y_font_size);
        else
            yticks([]);
        end
        hold off;
    end
end
end