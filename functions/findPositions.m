function positions = findPositions(data,Y,N)

if nargin < 3
    N = 1;
end

if nargin < 2
    Y = 1;
end
    % Initialize the matrix to store the positions for each row
    % Columns are: [first 25% 50% 75% last]
    positions = NaN(size(data, 1), 5);

    % Loop through each row (each animal)
    for i = 1:size(data, 1)
        % Find indices of Ns in the current row
        indicesOfN = find(data(i, :) == N);

        % Check if there are any Ns
        if ~isempty(indicesOfN)
            if length(indicesOfN) >= Y
            % Store the first occurrence
            positions(i, 1) = indicesOfN(1);                            % First
            % Calculate the positions for 25%, 50%, 75%, and last
            positions(i, 2) = indicesOfN(ceil(0.25 * length(indicesOfN))); % 25%
            positions(i, 3) = indicesOfN(ceil(0.50 * length(indicesOfN))); % 50%
            positions(i, 4) = indicesOfN(ceil(0.75 * length(indicesOfN))); % 75%
            positions(i, 5) = indicesOfN(end);                            % Last
            end
        end
    end
end
