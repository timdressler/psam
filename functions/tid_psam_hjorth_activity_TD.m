function activity = tid_psam_hjorth_activity_TD(signal)
% tid_psam_hjorth_activity_TD.m
% 
% Computes the Hjorth Activity.
%
% Usage:
%   tid_psam_hjorth_activity_TD(data)
%
% Inputs:
%   data - a numerical array [1-by-1], [Channels-by-1] or [Channels-by-Epochs]
%
% Outputs:
%   acitivity - the Hjorth Activity [1-by-1], [Channels-by-1] or [Channels-by-Epochs]
%
% The Hjorth Activity is defined as the variance of the signal.
%   For details, see Alawee et al. (2023). 
%   Note: In this implementation, variance is consistently calculated using 'n' as the denominator, 
%   rather than 'n - 1' as suggested by Alawee et al. (2023). 
%   Alawee et al. (2023) are inconsistent regarding whether 'n' or 'n - 1' 
%   should be used. Here, we adopt the use of 'n' throughout for consistency.
%
% Literature
    % Alawee, W. H., Basem, A., & Al-Haddad, L. A. (2023). 
    %   Advancing biomedical engineering: Leveraging Hjorth features for electroencephalography signal analysis. 
    %   Journal of Electrical Bioimpedance, 14(1), 66. https://doi.org/10.2478/joeb-2023-0009
%
% Tim Dressler, 06.05.2025

dims = ndims(signal);

switch dims
    case 2  % [1-by-Samples] or [Channels-by-Samples]
        activity = var(signal, 1, 2); % Result: [1-by-1] or [Channels-by-1]
    case 3  % [Channels-by-Samples-by-Epochs]
        activity = squeeze(var(signal, 1, 2));  % Result: [Channels-by-Epochs]
    otherwise
        error('Unsupported signal dimensions. Expected 2D or 3D array.');
end
end
