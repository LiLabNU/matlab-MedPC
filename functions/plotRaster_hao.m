function h = plotRaster_hao(data, type, color, s, ax)
% INPUT:
% data: A matrix where each row represents a trial and each column
%       represents a time point. A value of 1 represents an event,
%       and a value of 0 represents no event.

if nargin < 4
    color = 'k';
end

% Set the current axes to the provided ax for plotting
axes(ax);

hold on;
if type == "dot"
    [numTrials, ~] = size(data);
    h = cell(numTrials, 1);
    for trial = 1:numTrials
        spikeTimes = find(data(trial, :) == 1);
        h{trial} = plot(spikeTimes, trial * ones(size(spikeTimes)), color, 'MarkerSize', s);
    end
elseif type == "line"
    [numTrials, ~] = size(data);
    h = cell(numTrials, 1);
    for trial = 1:numTrials
        spikeTimes = find(data(trial, :) == 1);
        for t = spikeTimes
            h{trial} = plot([t, t], [trial-0.4, trial+0.4], color, 'LineWidth', s);
        end
    end
end
hold off;
end