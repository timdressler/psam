function transition_bandwidth = tid_psam_get_transition_bandwidth(cuttoff_freq)
% tid_psam_get_transition_bandwidth computes transition bandwidth based on the following rule.
%
% Assumes Hamming window!
%
% Transition bandwidth:
%   cutoff frequency <= 1 Hz: twice cutoff frequency
%   cutoff frequency > 1 Hz && <= 8 Hz: 2 Hz
%   cutoff frequency > 8: 25% of cutoff frequency
%
% Usage:
%   tid_psam_get_transition_band(cuttoff_freq)
%
% Inputs:
%   cuttoff_freq - cutoff frequency (e.g. 1 if 1 Hz is used)
%
% Outputs:
%   transition_bandwidth - calculated transition bandwidth
%
% Rule developed by Mareike Daeglau based on https://eeglab.org/others/Firfilt_FAQ.html
%
% Tim Dressler, 16.06.2024

if cuttoff_freq <= 1
    transition_bandwidth = 2 * cuttoff_freq;
elseif cuttoff_freq <= 8 && cuttoff_freq > 1
    transition_bandwidth = 2;
elseif cuttoff_freq > 8
    transition_bandwidth = 0.25 * cuttoff_freq;
end

disp(['Calculated transition band width: ' num2str(transition_bandwidth)])

end


