% bemobil_detect_bad_channels_modified_TD - Repeatedly detects bad channels using the "clean_artifacts" function of the EEGLAB
% "clean_rawdata" plugin. Computes average reference before detecting bad channels and uses an 0.5Hz highpass filter
% cutoff (inside clan_artifacts). Plots the rejections of each iteration and the final rejections. Plots segments of the
% data to check the rejected channels.
%
% Usage:
%   >>  [chans_to_interp, chan_detected_fraction_threshold, detected_bad_channels, rejected_chan_plot_handle, detection_plot_handle] =...
%           bemobil_detect_bad_channels(EEG, ALLEEG, CURRENTSET, chancorr_crit, chan_max_broken_time, chan_detect_num_iter,...
%           chan_detected_fraction_threshold, num_chan_rej_max_target, flatline_crit, line_noise_crit)
% 
% Inputs:
%   EEG                                 - current EEGLAB EEG structure
%   ALLEEG                              - complete EEGLAB data set structure
%   CURRENTSET                          - index of current EEGLAB EEG structure within ALLEEG
%   chancorr_crit                       - Correlation threshold. If a channel is correlated at less than this value
%                                           to its robust estimate (based on other channels), it is considered abnormal in
%                                           the given time window. OPTIONAL, default = 0.8.
%   chan_max_broken_time                - Maximum time (either in seconds or as fraction of the recording) during which a 
%                                           retained channel may be broken. Reasonable range: 0.1 (very aggressive) to 0.6
%                                           (very lax). OPTIONAL, default = 0.5.
%   chan_detect_num_iter                - Number of iterations the bad channel detection should run (default = 10)
%   chan_detected_fraction_threshold    - Fraction how often a channel has to be detected to be rejected in the final
%                                           rejection (default 0.5)
%   num_chan_rej_max_target             - Hard limit for the maximum amount of channel rejection. If the number of 
%                                           detected bad channels exceeds this limit, a 3-tier tie-breaker is applied to 
%                                           strictly select the worst channels: 1) Detection frequency across iterations 
%                                           ("badness"), 2) Flatline status (variance near zero), 3) Highest overall 
%                                           signal variance. If empty, use only chan_detected_fraction_threshold. 
%                                           Can be either a fraction of all channels (will be rounded, e.g. 1/5 of chans) 
%                                           or a specific integer number. (default = 1/5)
%   flatline_crit                       - Maximum duration a channel can be flat in seconds (default 'off')
%   line_noise_crit                     - If a channel has more line noise relative to its signal than this value, in
%                                           standard deviations based on the total channel population, it is considered
%                                           abnormal. (default: 'off')
%
% Outputs:
%   chans_to_interp                     - vector with channel indices to remove
%   chan_detected_fraction_threshold    - 'badness' threshold used for bad channel selection. Might deviate from 
%                                           chan_detected_fraction_threshold when there are a lot of bad channels.
%   detected_bad_channels               - n*m matrix of boolean indicating what channel was marked for rejection in which iteration
%                                           n is number of channels and m number of iterations
%   rejected_chan_plot_handle           - handle to the plot of the data segments to check cleaning
%   detection_plot_handle               - handle to the plot of the rejection iterations
%   
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%   EEGLAB, bemobil_avref, clean_artifacts
%
% Authors: Lukas Gehrke, 2017, Marius Klug, 2021, Timotheus Berg, 2021
%          Tim Dressler, 2026 (Modified 12.06.26: Enforced hard limit for num_chan_rej_max_target via 3-tier tie-breaker)

function [chans_to_interp, chan_detected_fraction_threshold, detected_bad_channels, rejected_chan_plot_handle, detection_plot_handle] = bemobil_detect_bad_channels_modified_TD(EEG, ALLEEG, CURRENTSET,...
    chancorr_crit, chan_max_broken_time, chan_detect_num_iter, chan_detected_fraction_threshold, num_chan_rej_max_target, flatline_crit, line_noise_crit)

if ~exist('chancorr_crit','var') || isempty(chancorr_crit)
	chancorr_crit = 0.8;
end
if ~exist('chan_max_broken_time','var') || isempty(chan_max_broken_time)
	chan_max_broken_time = 0.5;
end
if ~exist('chan_detect_num_iter','var') || isempty(chan_detect_num_iter)
	chan_detect_num_iter = 10;
end
if ~exist('chan_detected_fraction_threshold','var') || isempty(chan_detected_fraction_threshold)
	chan_detected_fraction_threshold = 0.5;
end
if ~exist('flatline_crit','var') || isempty(flatline_crit)
	flatline_crit = 'off';
end
if ~exist('line_noise_crit','var') || isempty(line_noise_crit)
	line_noise_crit = 'off';
end
if ~exist('num_chan_rej_max_target','var') || isempty(num_chan_rej_max_target)
	num_chan_rej_max_target = 1/5;
end

rejected_chan_plot_handle = [];
detection_plot_handle = [];

%%
if ~strcmp(EEG.ref,'average')
    disp('Re-referencing ONLY for bad channel detection now.')
    % compute average reference before finding bad channels 
    [ALLEEG, EEG, CURRENTSET] = bemobil_avref( EEG , ALLEEG, CURRENTSET);
end

%%
disp('Repeated bad channels detection...')
detected_bad_channels = [];
for i = 1:chan_detect_num_iter
    disp(['Iteration ' num2str(i) '/' num2str(chan_detect_num_iter)])
    clear hlp_microcache
    % remove bad channels, use default values of clean_artifacts, but specify just in case they may change
    [EEG_chan_removed,EEG_highpass,~,detected_bad_channels(1:EEG.nbchan,i)] = clean_artifacts(EEG,...
        'burst_crit','off','window_crit','off','channel_crit_maxbad_time',chan_max_broken_time,...
        'chancorr_crit',chancorr_crit,'line_crit',line_noise_crit,'highpass_band',[0.25 0.75],'flatline_crit',flatline_crit);

end
disp('...iterative bad channel detection done!')

badness_percent = sum(detected_bad_channels,2) / size(detected_bad_channels,2);

% Initial pass: identify bad channels based on the threshold
chans_to_interp = badness_percent >= chan_detected_fraction_threshold;

% Enforce the hard limit with a 3-tier tie-breaker
if ~isempty(num_chan_rej_max_target)
    if num_chan_rej_max_target < 1 && num_chan_rej_max_target > 0
        num_chan_rej_max_target = round(EEG.nbchan * num_chan_rej_max_target);
    end
    
    if sum(chans_to_interp) > num_chan_rej_max_target
        % --- Calculate Secondary Metrics ---
        % 1. Variance (Dim 2 is time)
        chan_variance = var(EEG.data, 0, 2); 
        
        % 2. Flatline check (Variance is practically zero)
        % Using 1e-10 to account for tiny floating-point inaccuracies
        is_flatline = double(chan_variance < 1e-10); 
        
        % Combine into a 3-column matrix: [Badness, Is_Flatline, Variance]
        sort_metrics = [badness_percent, is_flatline, chan_variance];
        
        % sortrows sorts by Col 1 first, then Col 2, then Col 3.
        % 'descend' ensures we get highest badness -> flatlines first -> highest noise.
        [~, sort_idx] = sortrows(sort_metrics, [1, 2, 3], 'descend');
        
        % Update the threshold variable for your output logs
        chan_detected_fraction_threshold = badness_percent(sort_idx(num_chan_rej_max_target));
        
        % Reset and apply hard limit based on the new intelligent sort
        chans_to_interp = false(size(badness_percent));
        chans_to_interp(sort_idx(1:num_chan_rej_max_target)) = true;
    end
end


%% select the final channels to remove and remove them from a dataset to plot

% give actual channel numbers as output
chans_to_interp = find(chans_to_interp);

disp('Detected bad channels: ')
disp({EEG.chanlocs(chans_to_interp).labels})

% take EOG out of the channels to interpolate: EOG is very likely to be different from the others but rightfully so
disp('Ignoring EOG channels for interpolation:')

disp({EEG.chanlocs(chans_to_interp(strcmpi({EEG.chanlocs(chans_to_interp).type},'EOG'))).labels})
chans_to_interp(strcmp({EEG.chanlocs(chans_to_interp).type},'EOG'))=[];

disp('Final bad channels: ')
disp({EEG.chanlocs(chans_to_interp).labels})


% remove channels and store channel mask for plotting
EEG_chan_removed = pop_select( EEG_highpass,'nochannel',chans_to_interp);

clean_channel_mask = ones(EEG.nbchan,1);
clean_channel_mask(chans_to_interp) = 0;
EEG_chan_removed.etc.clean_channel_mask = logical(clean_channel_mask);

