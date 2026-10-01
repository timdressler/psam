function tid_psam_save_params(pipeline_params, script_name, outpath)
% tid_psam_save_params.m
%
% Dynamically reads a parameter structure and saves it as a JSON 
% log file in the specified output folder.
%
% Usage:
%    tid_psam_save_params(pipeline_params, 'tid_psam_ica_preprocessing', OUTPATH)
%
% Inputs:
%    pipeline_params - Structure containing the parameters to log.
%    script_name     - String or char array of the script name (used for filename).
%    outpath         - String or char array of the directory to save the log.
%
% Outputs:
%    None (Writes a .json file to disk).
%
% Tim Dressler, 11.06.2026

% Check for valid inputs
if ~isstruct(pipeline_params)
    error('Input pipeline_params must be a structure.');
end

if ~ischar(script_name) && ~isstring(script_name)
    error('Input script_name must be a string or character array.');
end

if ~isfolder(outpath)
    error('Output path "%s" does not exist.', outpath);
end

% Define file ID and open file (note the .json extension)
log_filename = fullfile(outpath, sprintf('%s_params.json', char(script_name)));
fileID = fopen(log_filename, 'w');

if fileID == -1
    error('Could not open file for writing: %s', log_filename);
end

% Encode structure to a JSON-formatted string
% 'PrettyPrint' makes the JSON human-readable with line breaks and indents 
% (Requires MATLAB R2021a or newer. If using an older version, remove the argument.)
try
    json_string = jsonencode(pipeline_params, 'PrettyPrint', true);
catch
    % Fallback for older MATLAB versions that don't support 'PrettyPrint'
    json_string = jsonencode(pipeline_params);
end

% Write JSON string to file
fprintf(fileID, '%s', json_string);

% Close file and notify user
fclose(fileID);
disp(['Parameters successfully logged to: ', log_filename]);

end