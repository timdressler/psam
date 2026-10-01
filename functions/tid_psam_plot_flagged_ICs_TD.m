function tid_psam_plot_flagged_ICs_TD(EEG, sgtitlestr, varargin)
% tid_psam_plot_flagged_ICs_TD. 
%
% Plots extracted ICs line-by-line into a scrollable HTML report.
%
% Usage:
%   tid_psam_plot_flagged_ICs_TD(EEG, sgtitlestr, 'SavePath', savePath)
%   tid_psam_plot_flagged_ICs_TD(EEG, sgtitlestr, 'SavePath', savePath, 'PlotOn', false)
%
% Inputs:
%   EEG        - an EEGLAB structure containing information regarding flagged ICs 
%                and ICLabel classifications.
%   sgtitlestr - a string used for the main title of the HTML report.
%
% Required Name-Value Inputs:
%   'SavePath' - Full file path to save the HTML report. Must be a non-empty
%                string or character array (e.g., '/path/to/report.html').
%
% Optional Name-Value Inputs:
%   'PlotOn'        - Logical flag. If true, opens the generated HTML file in 
%                     your system's default web browser (default: true).
%   'KeepPlotAlive' - Maintained for backward compatibility (ignored).

% Parse optional inputs
p = inputParser;
p.FunctionName = 'tid_psam_plot_flagged_ICs_TD';
addRequired(p,  'EEG');
addRequired(p,  'sgtitlestr');
addParameter(p, 'SavePath',      '',   @(x) ischar(x) || isstring(x));
addParameter(p, 'PlotOn',        true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'KeepPlotAlive', true, @(x) islogical(x) && isscalar(x));
parse(p, EEG, sgtitlestr, varargin{:});

savePath = p.Results.SavePath;
if isempty(savePath)
    error('tid_psam_plot_flagged_ICs_TD: ''SavePath'' is mandatory and must be a non-empty string.');
end

plotOn = p.Results.PlotOn;

% Enforce .html extension for the path
[pathstr, name, ~] = fileparts(savePath);
if isempty(pathstr)
    savePath = fullfile(pwd, [name '.html']);
else
    savePath = fullfile(pathstr, [name '.html']);
end

% Check for flagged ICs data
if ~isfield(EEG.reject, 'gcompreject') || isempty(EEG.reject.gcompreject)
    error('No information regarding flagged ICs found in EEG.reject.gcompreject.');
end

% Check for ICLabel data
if ~isfield(EEG, 'etc') || ~isfield(EEG.etc, 'ic_classification') || ...
   ~isfield(EEG.etc.ic_classification, 'ICLabel')
    error('ICLabel classification data not found in EEG.etc.ic_classification.ICLabel.');
end

% Standard ICLabel categories
categories = {'Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise', 'Channel Noise', 'Other'};
ic_classes = EEG.etc.ic_classification.ICLabel.classifications;
numICs = size(EEG.icaweights, 1);

% Open HTML file for writing
fid = fopen(savePath, 'w');
if fid == -1
    error('Could not open file for writing: %s', savePath);
end

% Write HTML Header & CSS Styles
fprintf(fid, '<!DOCTYPE html>\n<html>\n<head>\n<title>%s</title>\n', sgtitlestr);
fprintf(fid, '<style>\n');
fprintf(fid, '  body { font-family: Arial, sans-serif; margin: 30px; background-color: #f7f9fa; color: #333; }\n');
fprintf(fid, '  h1 { text-align: center; color: #2c3e50; margin-bottom: 30px; }\n');
fprintf(fid, '  .ic-container { display: flex; align-items: center; background: white; border: 1px solid #e1e8ed; border-radius: 8px; margin-bottom: 20px; padding: 15px; box-shadow: 0 2px 4px rgba(0,0,0,0.04); }\n');
fprintf(fid, '  .ic-container.flagged { border-left: 8px solid #e74c3c; background-color: #fdf2f2; }\n');
fprintf(fid, '  .plot-col { flex: 0 0 250px; text-align: center; }\n');
fprintf(fid, '  .info-col { flex: 1; padding-left: 30px; }\n');
fprintf(fid, '  .ic-title { font-size: 1.3em; font-weight: bold; margin-bottom: 10px; color: #2c3e50; }\n');
fprintf(fid, '  .badge { display: inline-block; padding: 4px 8px; font-size: 0.85em; font-weight: bold; border-radius: 4px; text-transform: uppercase; }\n');
fprintf(fid, '  .badge.yes { background-color: #e74c3c; color: white; }\n');
fprintf(fid, '  .badge.no { background-color: #2ecc71; color: white; }\n');
fprintf(fid, '  .metrics-table { width: 100%%; margin-top: 10px; border-collapse: collapse; }\n');
fprintf(fid, '  .metrics-table td { padding: 4px 8px; font-size: 0.9em; }\n');
fprintf(fid, '  .bar-container { background-color: #e0e0e0; border-radius: 4px; width: 100%%; height: 12px; position: relative; }\n');
fprintf(fid, '  .bar { background-color: #3498db; height: 100%%; border-radius: 4px; }\n');
fprintf(fid, '  .bar.high { background-color: #e67e22; }\n');
fprintf(fid, '</style>\n</head>\n<body>\n');
fprintf(fid, '<h1>%s</h1>\n', sgtitlestr);

% Create an invisible temporary figure for plotting
hidden_fig = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 400 400]);

fprintf('Generating HTML report for %d ICs...\n', numICs);

% Loop through each independent component
for i = 1:numICs
    % 1. Generate the Topoplot invisibly
    clf(hidden_fig);
    set(0, 'CurrentFigure', hidden_fig); % <--- FORCES focus to this exact figure
    ax = axes('Parent', hidden_fig);     % <--- Creates explicit axes
    topoplot(EEG.icawinv(:, i), EEG.chanlocs, 'numcontour', 0, 'electrodes', 'off');
    set(ax, 'Position', [0.05 0.05 0.9 0.9]); % <--- Uses the explicit axis handle
    
    % 2. Save plot to a temporary image and convert to Base64 string
    tmp_img_path = [tempname '.png'];
    exportgraphics(hidden_fig, tmp_img_path, 'Resolution', 120);
    
    img_fid = fopen(tmp_img_path, 'r');
    img_bytes = fread(img_fid, Inf, '*uint8');
    fclose(img_fid);
    delete(tmp_img_path); 
    
    base64_str = char(matlab.net.base64encode(img_bytes));
    img_src = ['data:image/png;base64,', base64_str];
    
    % 3. Extract metadata info
    is_flagged = EEG.reject.gcompreject(i) == 1;
    if is_flagged
        row_class = 'ic-container flagged';
        badge_html = '<span class="badge yes">FLAGGED: YES</span>';
    else
        row_class = 'ic-container';
        badge_html = '<span class="badge no">FLAGGED: NO</span>';
    end
    
    % 4. Start writing the row container
    fprintf(fid, '<div class="%s">\n', row_class);
    
    % Left Column: The image
    fprintf(fid, '  <div class="plot-col">\n');
    fprintf(fid, '    <div class="ic-title">IC %d</div>\n', i);
    fprintf(fid, '    <img src="%s" width="220" height="220" alt="IC %d Topoplot">\n', img_src, i);
    fprintf(fid, '  </div>\n');
    
    % Right Column: Metadata & ICLabel Probabilities
    fprintf(fid, '  <div class="info-col">\n');
    fprintf(fid, '    <div style="margin-bottom: 12px;">%s</div>\n', badge_html);
    fprintf(fid, '    <table class="metrics-table">\n');
    
    % Loop through categories to populate progress bars
    for cat_idx = 1:7
        cat_name = categories{cat_idx};
        prob = ic_classes(i, cat_idx);
        prob_pct = prob * 100;
        
        bar_class = 'bar';
        if prob > 0.5
            bar_class = 'bar high';
        end
        
        fprintf(fid, '      <tr>\n');
        fprintf(fid, '        <td style="width: 120px; font-weight: bold;">%s:</td>\n', cat_name);
        fprintf(fid, '        <td style="width: 60px; text-align: right; padding-right: 15px;">%.1f%%</td>\n', prob_pct);
        fprintf(fid, '        <td>\n');
        fprintf(fid, '          <div class="bar-container">\n');
        fprintf(fid, '            <div class="%s" style="width: %.1f%%;"></div>\n', bar_class, prob_pct);
        fprintf(fid, '          </div>\n');
        fprintf(fid, '        </td>\n');
        fprintf(fid, '      </tr>\n');
    end
    
    fprintf(fid, '    </table>\n');
    fprintf(fid, '  </div>\n'); 
    fprintf(fid, '</div>\n\n'); 
end

% Close document structure
fprintf(fid, '</body>\n</html>\n');
fclose(fid);
close(hidden_fig);

fprintf('HTML Report successfully saved to: %s\n', savePath);

% 5. Auto-open report in web browser if requested
if plotOn
    web(savePath, '-browser');
end

end
