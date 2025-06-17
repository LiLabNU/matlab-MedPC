function [couples] = fillEmptytoSturcture(couples,emptyRowIndices)

% Get field names from the existing structure
fields = fieldnames(couples);

% Create an empty structure with the same fields as couples
emptyEntry = cell2struct(cell(size(fields)), fields, 1);

% Insert empty rows
for idx = numel(emptyRowIndices):-1:1
    insertIdx = emptyRowIndices(idx);
    couples = [couples(1:insertIdx-1); emptyEntry; couples(insertIdx:end)];
end
end