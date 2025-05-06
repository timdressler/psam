function activity = tid_psam_hjorth_activity_TD(signal)
% tid_psam_hjorth_activity_TD - Computes the Hjorth Activity.
%
% Usage:
%   tid_psam_hjorth_activity_TD(data)
%
% Inputs:
%   data - a 1-by-Samples, Channels-by-Samples or a Channels-by-Samples-by-Epochs numerial array
%
% Outputs:
%   acitivity - the Hjorth Activity [1-by-1], [Channels-by-1] or [Channels-by-Epochs]
%
% The Hjorth Activity is defined as the variance of the signal.
%
% Tim Dressler, 07.11.2024

dims = ndims(signal);

switch dims
    case 2  % [1-by-Samples] or [Channels-by-Samples]
        activity = var(signal, [], 2); % Result: [1-by-1] or [Channels-by-1]
    case 3  % [Channels-by-Samples-by-Epochs]
        activity = squeeze(var(signal, 0, 2));  % Result: [Channels-by-Epochs]
    otherwise
        error('Unsupported signal dimensions. Expected 2D or 3D array.');
end
end
