function [averages] = segmentTimestamps(referenceTimestamps, applyingTimestamps, values, numSegments, segmentType)
    if ~isnumeric(referenceTimestamps) || ~isnumeric(applyingTimestamps) || ~isnumeric(values)
        error('All input data must be numeric arrays.');
    end

    if strcmp(segmentType, 'position')
        % Position-based segmentation logic (existing logic)
        referenceTimestamps = sort(referenceTimestamps);
        minTime = referenceTimestamps(1);
        maxTime = referenceTimestamps(end);
        totalRange = maxTime - minTime;
        segmentSize = totalRange / (numSegments - 1);

        segments = [minTime + (0:numSegments-2)' * segmentSize, minTime + (1:numSegments-1)' * segmentSize - 1; ...
                    minTime + (numSegments-1) * segmentSize, maxTime];
    elseif strcmp(segmentType, 'event')
        % Event-based segmentation logic
        if length(referenceTimestamps) < numSegments
            error('Not enough events to segment.');
        end
        referenceTimestamps = sort(referenceTimestamps);
        eventsPerSegment = round(length(referenceTimestamps) / numSegments);
        segments = zeros(numSegments, 2);
        for i = 1:numSegments
            if i == 1
                segments(i, 1) = referenceTimestamps(1);
            else
                if (i-1) * eventsPerSegment + 1 <= length(referenceTimestamps)
                    segments(i, 1) = referenceTimestamps((i-1) * eventsPerSegment + 1);
                else
                    segments(i, 1) = referenceTimestamps(end);
                end
            end
            
            if i == numSegments
                segments(i, 2) = referenceTimestamps(end);
            else
                segments(i, 2) = referenceTimestamps(i * eventsPerSegment);
            end
        end
    else
        error('Invalid segment type specified. Use ''position'' or ''event''.');
    end

    lowerBounds = repmat(segments(:, 1)', length(applyingTimestamps), 1);
    upperBounds = repmat(segments(:, 2)', length(applyingTimestamps), 1);
    timestamps = repmat(applyingTimestamps, 1, numSegments);

    inSegment = timestamps >= lowerBounds & timestamps <= upperBounds;
    [~, ids] = max(inSegment, [], 2);

    sums = accumarray(ids, values, [numSegments, 1]);
    counts = accumarray(ids, 1, [numSegments, 1]);

    averages = zeros(numSegments, 1);
    valid = counts > 0;
    averages(valid) = sums(valid) ./ counts(valid);
end
