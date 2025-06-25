function mobility = tid_psam_hjorth_mobility_TD(signal)
% tid_psam_hjorth_mobility_TD - Computes the Hjorth Mobility.
%
% Usage:
%   tid_psam_hjorth_mobility_TD(data)
%
% Inputs:
%   data - a numerical array [1-by-1], [Channels-by-1] or [Channels-by-Epochs]
%
% Outputs:
%   mobility - the Hjorth Mobility [1-by-1], [Channels-by-1] or [Channels-by-Epochs]
%
% The Hjorth Mobility is defined as the square root of variance of first derivative divided by variance of the signal.
%   For details, see Alawee et al. (2023). 
%   Note: In this implementation, variance is consistently calculated using 'n' as the denominator, 
%   rather than 'n - 1' as suggested by Alawee et al. (2023). 
%   Alawee et al. (2023) are inconsistent regarding whether 'n' or 'n - 1' 
%   should be used. Here, we adopt the use of 'n' throughout for consistency.
%
% Also see: tid_psam_hjorth_activity_TD.
%
% Literature
    % Alawee, W. H., Basem, A., & Al-Haddad, L. A. (2023). 
    %   Advancing biomedical engineering: Leveraging Hjorth features for electroencephalography signal analysis. 
    %   Journal of Electrical Bioimpedance, 14(1), 66. https://doi.org/10.2478/joeb-2023-0009
%
% Tim Dressler, 06.05.2025

    diff1 = diff(signal, 1, 2);

    var0 = tid_psam_hjorth_activity_TD(signal);
    var1 = tid_psam_hjorth_activity_TD(diff1);

    mobility = sqrt(var1 ./ var0);
end
