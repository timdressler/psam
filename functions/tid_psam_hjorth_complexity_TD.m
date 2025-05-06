function complexity = tid_psam_hjorth_complexity_TD(signal)
% tid_psam_hjorth_complexity_TD - Computes the Hjorth Complexicty.
%
% Usage:
%   tid_psam_hjorth_complexity_TD(data)
%
% Inputs:
%   data - a numerical array [1-by-1], [Channels-by-1] or [Channels-by-Epochs]
%
% Outputs:
%   complexity - the Hjorth Complexity [1-by-1], [Channels-by-1] or [Channels-by-Epochs]
%
% The Hjorth Mobility is defined as the square root of variance of first derivative divided by variance of the signal.
%   For details, see Alawee et al. (2023)
%
% Also see: tid_psam_hjorth_activity_TD, tid_psam_hjorth_mobility_TD.
%
% Literature
    % Alawee, W. H., Basem, A., & Al-Haddad, L. A. (2023). 
    %   Advancing biomedical engineering: Leveraging Hjorth features for electroencephalography signal analysis. 
    %   Journal of Electrical Bioimpedance, 14(1), 66. https://doi.org/10.2478/joeb-2023-0009
%
% Tim Dressler, 06.05.2025

    diff1 = diff(signal, 1, 2);
    
    mobility = tid_psam_hjorth_mobility_TD(signal);
    mobility_deriv = tid_psam_hjorth_mobility_TD(diff1);

    complexity = mobility_deriv ./ mobility;
end
