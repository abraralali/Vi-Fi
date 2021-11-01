function NewzedBoxgndmatch = importHansiMatchfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  NEWZEDBOXGNDMATCH = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  NEWZEDBOXGNDMATCH = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  NewzedBoxgndmatch = importfile("C:\RANProject\WINLABData\20211007_113104\NewzedBox_gnd_match.txt", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 24-Oct-2021 17:34:37

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 11);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
% opts.VariableNames = ["png", "VarName2", "Nicholas", "VarName4", "Bo", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11"];
opts.VariableTypes = ["string", "double", "categorical", "double", "categorical", "double", "categorical", "double", "categorical", "double", "categorical"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";

% Specify variable properties
% opts = setvaropts(opts, "png", "WhitespaceRule", "preserve");
% opts = setvaropts(opts, ["png", "Nicholas", "Bo", "VarName7", "VarName9", "VarName10", "VarName11"], "EmptyFieldRule", "auto");

% Import the data
NewzedBoxgndmatch = readtable(filename, opts);

end