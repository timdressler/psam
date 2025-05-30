function activity = tid_psam_plot_flagged_ICs_TD(EEG,sgtitlestr,savestr)
% tid_psam_plot_flagged_ICs_TD - Plots flagged ICs..
%
% Usage:
%   tid_psam_plot_flagged_ICs_TD(EEG, sgtitlestr, savestr)
%
% Inputs:
%   EEG - an EEGLAB structure
%       EEG must contain information regarding flagged ICs
%
%   sgtitlestr - a string used for the title of the plot
%
%   savestr - a sting used for saving the plot (full path required)
%
% Outputs:
%   -
%
% Plots:
%   Plots all extracted IC and tags the flagged ones.
%
% Tim Dressler, 24.05.2025

% Check if ICLabel Plugin (Pion-Tonachini et al., 2019) was applied and the ICs are flagged
if  ~isfield(EEG.reject,'gcompreject') || isempty(EEG.reject.gcompreject)
    error('No information regarding flagged ICs found.')
end

% Get flagged ICs
flagged_comps = find(EEG.reject.gcompreject);

% Create plot
figure('Units', 'normalized', 'Position', [0 0.0500 1 0.8771]);
for i = 1:size(EEG.icaweights,1)
    subplot(ceil(sqrt(size(EEG.icaweights,1))), ceil(sqrt(size(EEG.icaweights,1))), i);
    topoplot(EEG.icawinv(:,i), EEG.chanlocs, 'numcontour', 0);
    titlestr = sprintf('IC %d', i);
    if ismember(i, flagged_comps) % Mark flagged ICs as such
        titlestr = [titlestr, ' FLAGGED'];  
    end
    title(titlestr);
end
sgtitle(sgtitlestr)

% Save figure
saveas(gcf,savestr)

disp('Plot saved!')

