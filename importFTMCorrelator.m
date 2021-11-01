function distances = importFTMCorrelator(filename, dataLines)
%IMPORTFILE Import data from a text file
%  DISTANCES = IMPORTFILE(FILENAME) reads data from text file FILENAME
%  for the default selection.  Returns the data as a table.
%
%  DISTANCES = IMPORTFILE(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  distances = importfile("C:\Users\aalali\CameraIMUCorrelation\Trials\NewTrials\distancesDTW\20201223_140951_distances.txt", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 22-Mar-2021 12:41:56

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4"];
opts.VariableTypes = ["datetime", "double", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "VarName4", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "VarName4", "EmptyFieldRule", "auto");
opts = setvaropts(opts, "VarName1", "InputFormat", "yyyy-MM-dd HH:mm:ss.SSS");

% Import the data
distances = readtable(filename, opts);

end