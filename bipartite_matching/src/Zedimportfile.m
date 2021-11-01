function zedBox = Zedimportfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  ZEDBOX = IMPORTFILE(FILENAME) reads data from text file FILENAME for
%  the default selection.  Returns the data as a table.
%
%  ZEDBOX = IMPORTFILE(FILE, DATALINES) reads data for the specified row
%  interval(s) of text file FILENAME. Specify DATALINES as a positive
%  scalar integer or a N-by-2 array of positive scalar integers for
%  dis-contiguous row intervals.
%
%  Example:
%  zedBox = importfile("C:\Users\aalali\CameraIMUCorrelation\Trials\NewTrials\20201223_140951\zedBox_20201223_140951.txt", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 28-Feb-2021 23:20:58

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["frameNo", "trackID", "bbx", "bby", "bbw", "bbh"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
zedBox = readtable(filename, opts);

end