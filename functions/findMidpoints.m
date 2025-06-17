function midpointPositions = findMidpoints(data,N)

if nargin < 2
    N = 1;
end
    % Initialize the array to store the midpoint positions for each row
    midpointPositions = zeros(size(data, 1), 1);

    % Loop through each row (each animal)
    for i = 1:size(data, 1)
        % Find indices of 1s in the current row
        indicesOfOnes = find(data(i, :) == N);

        % Calculate the midpoint index of the 1s if any exist
        if ~isempty(indicesOfOnes)
            midpointIndex = ceil(length(indicesOfOnes) / 2);
            % Get the actual position in the row
            midpointPositions(i) = indicesOfOnes(midpointIndex);
        else
            % If no 1s are present, set the position to NaN or a specific marker
            midpointPositions(i) = NaN;
        end
    end
end
