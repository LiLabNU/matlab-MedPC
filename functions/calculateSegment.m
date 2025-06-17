function sum_struct = calculateSegment(data, RefSeries, numSegments, position, segmentType)
% CALCULATESEGMENT Computes summed values over segments for each field in a data structure.
%
% INPUTS:
%   data         - A struct with fields that each contain a 2D matrix (trials × timepoints).
%   RefSeries    - A binary vector (1 × timepoints) used to define segmentation boundaries.
%                  If segmentType is 'time', 1s define the range for segmentation.
%                  If segmentType is 'event', 1s represent event timestamps for segment splitting.
%   numSegments  - Number of segments to divide the data into.
%   position     - Integer specifying the row of data.(fieldName) to analyze.
%   segmentType  - String, either 'time' or 'event'.
%                  'time' = split evenly up to the last event.
%                  'event' = split based on event occurrence counts.
%
% OUTPUT:
%   sum_struct   - A struct with the same fields as `data`. Each field contains a
%                  1×numSegments vector representing summed values for each segment.

% Initialize a structure to hold the mean values for each field
sum_struct = struct();

% Retrieve the names of all fields in the data structure
fieldNames = fieldnames(data);

% Loop through each field
for i = 1:length(fieldNames)
    fieldName = fieldNames{i};  % Get the field name
    DataSeries = data.(fieldName)(position,:);  % Access the field data

    % Use the updated function to calculate the segment sums
    sum_matrix = calculateSegmentSum(DataSeries, RefSeries, numSegments, segmentType);

    % Store the result in the corresponding field of the output structure
    sum_struct.(fieldName) = sum_matrix;
end
end


function sum_matrix = calculateSegmentSum(DataSeries, RefSeries, numSegments, segmentType)
if strcmp(segmentType, 'time')
    % Original position-based segmentation logic
    idx_last_one = find(RefSeries == 1, 1, 'last');
    if length(DataSeries) < idx_last_one
        error('The length of the Data series is less than the index of the last ''1'' in the Ref series.');
    end
    segment_length = floor(idx_last_one / numSegments);
    sum_matrix = zeros(1, numSegments);
    for i = 1:numSegments
        if i == numSegments
            segment = DataSeries((i-1)*segment_length + 1:idx_last_one);
        else
            segment = DataSeries((i-1)*segment_length + 1:i*segment_length);
        end
        sum_matrix(i) = sum(segment);
    end
elseif strcmp(segmentType, 'event')
    % Corrected event-based segmentation logic
    eventIndices = find(RefSeries == 1);
    if length(eventIndices) < numSegments
        error('Not enough events to segment.');
    end

    eventsPerSegment = round(length(eventIndices) / numSegments);
    sum_matrix = zeros(1, numSegments);
    for i = 1:numSegments

        if i == 1
            startIdx = 1;  % Start from the first data point for the first segment
        else
            startIdx = eventIndices((i-1) * eventsPerSegment) + 1;
        end

        if i == numSegments
            endIdx = length(DataSeries);  % Include all remaining data in the last segment
        else
            if (i * eventsPerSegment) > length(eventIndices)
                endIdx = length(DataSeries);  % Ensure not to exceed data length
            else
                endIdx = eventIndices(i * eventsPerSegment);
            end
        end

        segment = DataSeries(startIdx:endIdx);
        sum_matrix(i) = sum(segment);
    end
else
    error('Invalid segment type specified. Use ''time'' or ''event''.');
end
end