function [segmentedSum, miceInclude] = segmentEdges(refData, minSuc, minSegNum, currdata, currtrialTS, currFRdetails, currPEdetails, numSegTrials, mice)
% segmentedSumCalculation calculates segmentation summary statistics for each file.
%
%   Inputs:
%       refData      - Matrix (numFiles x timepoints) for reference events.
%       minSuc       - Minimum number of sucrose events required.
%       minEvent     - Minimum number of reference events required.
%       currdata     - Structure containing fields (e.g., Sucrose) for each file.
%       MedPCBin     - Binning factor to scale timestamps.
%       numSegments  - Number of segments to divide the time series.
%       currtrialTS  - Structure of time series (per event type) for each file.
%       currFRdetails- Cell array with FR details for each file.
%       currPEdetails- Cell array with PE details for each file.
%       refEvent     - String specifying the reference event name.
%       mice         - Array of mouse IDs.
%
%   Outputs:
%       segmentedSum - Structure with segmentation summary (each field is a numFiles x numSegments matrix).
%       miceInclude  - String array listing mouse IDs that passed inclusion criteria.

%% Initialization
numFiles = size(refData, 1);
segmentedSum = struct();
miceInclude = strings(numFiles, 1);

%% Process each file
for i = 1:numFiles
    RefSeries = refData(i, :);
    segNum = round(sum(RefSeries)/numSegTrials);
    includeMouse = segNum >= minSegNum && sum(currdata.Sucrose(i, :)) >= minSuc;

    if includeMouse
        miceInclude(i) = mice(i);
        edges = 0:numSegTrials:sum(RefSeries);
        edges = edges(2:end); % boundaries
        edgesT = find(RefSeries==1);
        edgesT = edgesT(edges);
        temp = {};
        for j = 1:segNum
            if j==1
                temp{j} = 1:edgesT(j);
            elseif j == segNum
                temp{j} = edgesT(j-1)+1:length(RefSeries);
            else
                temp{j} = edgesT(j-1)+1:edgesT(j);
            end
        end
        edgesTimestamps{i,:} = temp;
    end
end 

%% 
fieldNames = fieldnames(currdata);
numFields = length(fieldNames);
for k = 1:numFields+1
    if k <= numFields
        fieldName = fieldNames{k};
    else
        fieldName = "SegDuration";        
    end
    for i = 1:numFiles
        if numel(edgesTimestamps) >=i && ~isempty(edgesTimestamps{i})
            for j = 1:numel(edgesTimestamps{i})
                if k <= numFields
                    temp = currdata.(fieldName)(i,:);
                    tempSeg = temp(edgesTimestamps{i}{j});
                    if fieldName == "PortEntry"
                        segmentedSum.(fieldName)(i,j) = sum(tempSeg)/100;
                    else
                        segmentedSum.(fieldName)(i,j) = sum(tempSeg);
                    end
                else
                    segmentedSum.(fieldName)(i,j) = length(edgesTimestamps{i}{j})/100;
                end
            end
        else
            segmentedSum.(fieldName)(i,:) = zeros(1,size(segmentedSum.(fieldName),2));
        end
    end
    
end

fieldNames = fieldnames(segmentedSum);
for k = 1:numFields+1
    fieldName = fieldNames{k};
    for i = 1:numFiles
        n = find(segmentedSum.(fieldName)(i,:)==0);
        if ~isempty(n)
            segmentedSum.(fieldName)(i,n) = NaN;
        end
    end
end
%%



