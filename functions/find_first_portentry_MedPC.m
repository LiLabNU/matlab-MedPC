function firstPortEntries = find_first_portentry_MedPC(rewardTS, portEntryTS)
% Finds the first port entry after each reward based on reward positions.
%
% Inputs:
% - rewardTS: 1 × N numeric array of reward positions (indices)
% - portEntryTS: 1 × M binary array (1 = port entry, 0 = no port entry)
%
% Output:
% - firstPortEntries: 1 × N array containing the indices of the first port entry after each reward
%                     (NaN if no port entry is found after a reward)
%
% Example:
% rewardTS = [3, 7, 11]; % Positions of rewards
% portEntryTS = [0 1 0 0 0 1 0 0 1 0 0 1 0]; % Binary array of port entries
% firstPortEntries = find_first_portentry_bin_pos(rewardTS, portEntryTS);

    % Find indices of all port entries
    portEntryIdx = find(portEntryTS == 1);

    % Preallocate output
    firstPortEntries = nan(size(rewardTS));

    for i = 1:length(rewardTS)
        % Find the first port entry occurring after the current reward
        idx = find(portEntryIdx > rewardTS(i), 1, 'first');

        % Store the first port entry index (or NaN if none found)
        if ~isempty(idx)
            firstPortEntries(i) = portEntryIdx(idx);
        end
    end
end
