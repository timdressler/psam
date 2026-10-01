function tid_psam_unpack_params(pipeline_params)
% tid_psam_unpack_params.m
%
% Extracts all fields from a structure and assigns them as standalone 
% variables in the caller's workspace.
%
% Usage:
%    tid_psam_unpack_params(pipeline_params)
%
% Inputs:
%    pipeline_params - A structure containing parameter names and values.
%
% Outputs:
%    None (Variables are injected directly into the calling workspace).
%
% Tim Dressler, 11.06.2026

% Check for valid inputs
if ~isstruct(pipeline_params)
    error('Input pipeline_params must be a structure.');
end

% Handle empty struct case
if isempty(fieldnames(pipeline_params))
    warning('pipeline_params structure is empty. No variables unpacked.');
    return;
end

% Get all field names
fields = fieldnames(pipeline_params);

% Loop through each field and assign it to the caller workspace
for i = 1:length(fields)
    fname = fields{i};
    fval = pipeline_params.(fname);
    assignin('caller', fname, fval);
end

end