%% This function writes a given text in a log file as well as displays it
% on the matlab terminal
%
% Inputs:
% - fid points towards the file where to write
% - text is a string with the text to write
function [] = WriteInformation(fid,text)

    if ~isempty(fid)
        % Writes on log file (and goes to the next line for the following
        % command)
        fprintf(fid,text);
        fprintf(fid,'\n');
    end
    
    % Displays on terminal
    fprintf(text);
    fprintf('\n');
end