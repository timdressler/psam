function activity = tid_psam_hjorth_activity_TD(signal)
% tid_psam_hjorth_activity_TD - Computes the Hjorth Activity.
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
%   For details, see Alawee et al. (2023)
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
        activity = squeeze(var(signal, 0, 2));  % Result: [Channels-by-Epochs]
    otherwise
        error('Unsupported signal dimensions. Expected 2D or 3D array.');
end
end
