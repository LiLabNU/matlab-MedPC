function [segmentedSum, miceInclude] = segmentedSumCalculation(refData, minSuc, minEvent, currdata, MedPCBin, numSegments, currtrialTS, currFRdetails, currPEdetails, refEvent, mice)
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
fieldNames = fieldnames(currdata);
numFields = length(fieldNames);
numFiles = size(refData, 1);

segmentedSum = struct();
miceInclude = strings(numFiles, 1);

% Preallocate segmentation matrices for each field in currdata
for k = 1:numFields
    segmentedSum.(fieldNames{k}) = zeros(numFiles, numSegments);
end
% Additional fields
segmentedSum.interSucroseITI = zeros(numFiles, numSegments);
segmentedSum.numberOfInterSucrosePEs = zeros(numFiles, numSegments);
segmentedSum.interSucrosePEDuration = zeros(numFiles, numSegments);
segmentedSum.firstPELatency = zeros(numFiles, numSegments);
segmentedSum.firstPEDuration = zeros(numFiles, numSegments);
segmentedSum.FRCompletion = zeros(numFiles, numSegments);
segmentedSum.FRDuration = zeros(numFiles, numSegments);
segmentedSum.NPoverSuc = zeros(numFiles, numSegments);

%% Process each file
for i = 1:numFiles
    RefSeries = refData(i, :);
    includeMouse = sum(RefSeries) >= minEvent && sum(currdata.Sucrose(i, :)) >= minSuc;

    if includeMouse
        miceInclude(i) = mice(i);
        tempResults = calculateSegment(currdata, RefSeries, numSegments, i, 'event');
        for k = 1:numFields
            fieldName = fieldNames{k};
            segmentedSum.(fieldName)(i, :) = tempResults.(fieldName);
        end
    end
    
    % reference timestamps
    referenceTimeSeries = currtrialTS.(refEvent){i};

    % Calculate FR and sucrose details if available
    if ~isempty(currFRdetails{i})       
        segTimeSeries = currtrialTS.Sucrose{i};
        if includeMouse && ~isempty(segTimeSeries)
            % Calculate interSucroseITI, numberOfInterSucrosePEs, and interSucrosePEDuration
            dataSeries = currPEdetails{i}.interSucroseITI;
            if ~isempty(dataSeries)
                segmentedSum.interSucroseITI(i, :) = segmentTimestamps(referenceTimeSeries, segTimeSeries, dataSeries, numSegments, 'event');
            end
            dataSeries = currPEdetails{i}.numberOfInterSucrosePEs;
            if ~isempty(dataSeries)
                segmentedSum.numberOfInterSucrosePEs(i, :) = segmentTimestamps(referenceTimeSeries, segTimeSeries, dataSeries, numSegments, 'event');
            end
            dataSeries = currPEdetails{i}.interSucrosePEDuration;
            if ~isempty(dataSeries)
                segmentedSum.interSucrosePEDuration(i, :) = segmentTimestamps(referenceTimeSeries, segTimeSeries, dataSeries, numSegments, 'event');
            end
            % Calculate FRDuration
            dataSeries = currFRdetails{i}.FRDuration;
            if ~isempty(dataSeries)
                segmentedSum.FRDuration(i, :) = segmentTimestamps(referenceTimeSeries, segTimeSeries, dataSeries, numSegments, 'event');
            end
            % Calculate NP over Suc ratio
            dataSeries = currFRdetails{i}.nosepokesBetweenRewards;
            if ~isempty(dataSeries)
                segmentedSum.NPoverSuc(i, :) = segmentTimestamps(referenceTimeSeries, segTimeSeries, dataSeries, numSegments, 'event');
            end
        end
    end

%     if ~isempty(currPEdetails{i}) && ~isempty(currtrialTS.Shock{i})
%         referenceTimeSeries = currtrialTS.(refEvent){i};
%         segTimeSeries = currtrialTS.Shock{i};
%         % Calculate LatencyBetweenShockToPE
%         dataSeries = currPEdetails{i}.LatencyBetweenShockToPE;
%         if ~isempty(dataSeries) && ~isnan(dataSeries)
%             segmentedSum.LatencyBetweenShockToPE(i, :) = segmentTimestamps(referenceTimeSeries, segTimeSeries, dataSeries, numSegments, 'event');
%         end
%     end

    if ~isempty(currFRdetails{i})
        segTimeSeries = currPEdetails{i}.firstPEOnset;
        if includeMouse && ~isempty(segTimeSeries)
            % Calculate firstPELatency
            dataSeries = currPEdetails{i}.firstPELatency;
            if ~isempty(dataSeries)
                segmentedSum.firstPELatency(i, :) = segmentTimestamps(referenceTimeSeries, segTimeSeries, dataSeries, numSegments, 'event');
            end
            % Calculate firstPEDuration
            dataSeries = currPEdetails{i}.firstPEDuration;
            if includeMouse && ~isempty(currPEdetails{i}.eventTime) && ~isempty(dataSeries)
                segmentedSum.firstPEDuration(i, :) = segmentTimestamps(referenceTimeSeries, segTimeSeries, dataSeries, numSegments, 'event');
            end
        end
    end

    if ~isempty(currFRdetails{i})
        segTimeSeries = currtrialTS.FirstActive{i};
        if includeMouse && ~isempty(segTimeSeries)
            % Calculate FRCompletion
            dataSeries = currFRdetails{i}.FRcompletion;
            if ~isempty(dataSeries)
                segmentedSum.FRCompletion(i, :) = segmentTimestamps(referenceTimeSeries, segTimeSeries, dataSeries, numSegments, 'event');
            end
        end
    end
end

%% Final Adjustments: Scale results by MedPCBin and round up
segmentedSum.interSucroseITI = ceil(segmentedSum.interSucroseITI / MedPCBin);
segmentedSum.interSucrosePEDuration = ceil(segmentedSum.interSucrosePEDuration / MedPCBin);
segmentedSum.firstPELatency = ceil(segmentedSum.firstPELatency / MedPCBin);
segmentedSum.firstPEDuration = ceil(segmentedSum.firstPEDuration / MedPCBin);
segmentedSum.FRDuration = ceil(segmentedSum.FRDuration / MedPCBin);

% Note: segmentedSum.PortEntry is referenced below but must exist in currdata.
if isfield(segmentedSum, 'PortEntry')
    segmentedSum.PortEntry = ceil(segmentedSum.PortEntry / MedPCBin);
end
end
