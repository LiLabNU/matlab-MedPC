function [rBouts, uBouts, rewardComsumpRatio] = find_bouts(reward, portEntry, portExit, boutMaxGap, boutMinDur, currFolder, ch)
% FIND_BOUTS Identifies rewarded and unrewarded port entry bouts
%
% Inputs:
%   reward     - Reward timestamps
%   portEntry  - Port entry timestamps
%   portExit   - Port exit timestamps
%   boutMaxGap - Maximum gap between bouts to merge (in 10ms)
%   boutMinDur - Minimum bout duration to keep (in 10ms)
%   currFolder - Current folder name (for warning messages)
%   ch         - Channel number (for warning messages)
%
% Outputs:
%   rBouts - Rewarded bouts [start, end, duration]
%   uBouts - Unrewarded bouts [start, end, duration]

% Initialize output arrays
rBouts = [];
uBouts = [];

boutMaxGap = boutMaxGap*100;
boutMinDur = boutMinDur*100;

% Check for empty inputs
if isempty(portEntry) || isempty(portExit)
    warning('Empty portEntry or portExit timestamps for %s CH%d', currFolder, ch);
    return;
end

% If lengths don't match, take the minimum
minLength = min(length(portEntry), length(portExit));
if minLength < 1
    warning('Not enough entry/exit pairs for %s CH%d', currFolder, ch);
    return;
end

portEntry = portEntry(1:minLength);
portExit = portExit(1:minLength);

% Build initial bouts [start, end, duration]
bouts = [portEntry, portExit, portExit - portEntry];

% Filter out bouts that are too short
validIdx = bouts(:,3) >= boutMinDur;
bouts = bouts(validIdx, :);

% If no valid bouts after duration filtering, return empty arrays
if isempty(bouts)
    warning('No bouts meet the minimum duration criteria for %s CH%d', currFolder, ch);
    return;
end

% Merge bouts that are close together
if size(bouts, 1) > 1
    mergedBouts = bouts(1, :);  % Start with the first bout
    j = 1;  % Index for merged bouts

    for i = 2:size(bouts, 1)
        if (bouts(i, 1) - mergedBouts(j, 2)) <= boutMaxGap
            % Merge: extend end time and add actual duration
            mergedBouts(j, 2) = bouts(i, 2);  % New end time
            mergedBouts(j, 3) = mergedBouts(j, 3) + bouts(i, 3);  % Sum durations only
        else
            % Start a new bout
            j = j + 1;
            mergedBouts(j, :) = bouts(i, :);
        end
    end

    bouts = mergedBouts;
end


% Identify rewarded and unrewarded bouts
if isempty(reward)
    % All bouts are unrewarded if no rewards
    uBouts = bouts;
    return;
end

% For reward, check if there's a bout followed by reward delivery
isRewarded = false(size(bouts, 1), 1);
for r = 1:length(reward)
    % Find the first bout that starts after this reward
    idx = find(bouts(:,1) > reward(r), 1, 'first');
    if ~isempty(idx)
        isRewarded(idx) = true;
    end
end
rewardComsumpRatio = sum(isRewarded) / length(reward);

% Separate rewarded and unrewarded bouts
firstBoutsAfterRewards = [];  % To store valid bout onsets
boutOnsets = bouts(:,1);
for r = 1:length(reward)
    idx = find(boutOnsets > reward(r), 1, 'first');
    if ~isempty(idx)
        firstBoutsAfterRewards(end+1) = boutOnsets(idx);  % Append to result
    end
end

% Get unique bouts that were the first after a reward
rBouts = unique(firstBoutsAfterRewards(:));
uBouts = boutOnsets(~ismember(boutOnsets, rBouts));

end