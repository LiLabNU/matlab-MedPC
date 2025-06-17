function [N_control, N_exp] = plotTimeSeriesWithSEM_Multi(controlData, expData, enableStats)
% Function to plot control and experimental data across time points with SEM error bars
% and optionally perform repeated-measures two-way ANOVA and Bonferroni-corrected posthoc tests.
%
% Inputs:
%   - controlData: Matrix (m x n), m = mice, n = time points
%   - expData:     Matrix (m x n), same dimensions as controlData
%   - enableStats: true/false flag to perform stats and annotate plot

if nargin < 3
    enableStats = false;
end

% Sanity check
if size(controlData, 2) ~= size(expData, 2)
    error('controlData and expData must have the same number of columns (time points).');
end


% Remove mice (rows) that are all zeros
controlData = controlData(~all(controlData == 0, 2), :);
expData     = expData(~all(expData == 0, 2), :);

% Number of time points
numTimePoints = size(controlData, 2);

% Compute mean and SEM
meanControl = mean(controlData, 1, 'omitnan');
semControl  = std(controlData, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(controlData), 1));
meanExp     = mean(expData, 1, 'omitnan');
semExp      = std(expData, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(expData), 1));


% Get N for each group
N_control = size(controlData, 1);
N_exp     = size(expData, 1);
xVals = 1:numTimePoints;


% Plot
hold on;
errorbar(xVals, meanControl, semControl, '-o', 'LineWidth', 1.5, 'Color', 'b', 'DisplayName', 'Control');
errorbar(xVals, meanExp, semExp, '-o', 'LineWidth', 1.5, 'Color', 'r', 'DisplayName', 'Experimental');
ylabel('Mean ± SEM');
xticks(xVals);
grid on;

% Combine all Y values for bounds
allY = [meanControl + semControl, meanExp + semExp, meanControl - semControl, meanExp - semExp];
yMin = min(allY);
yMax = max(allY);
yRange = range(allY);

yLower = max(0, yMin - 0.05 * yRange);
yUpper = yMax + 0.25 * yRange;
ylim([yLower, yUpper]);

% Stats section
if enableStats
    

    % Remove rows with NaNs
    controlData = controlData(all(~isnan(controlData), 2), :);
    expData     = expData(all(~isnan(expData), 2), :);

    % Run RM mixed ANOVA using anova_rm
    try
        [p, ~] = anova_rm({controlData, expData}, 'off');
    catch ME
        warning('anova_rm failed: %s', char(ME.message));
        fprintf('Size of controlData: %s\n', mat2str(size(controlData)));
        fprintf('Size of expData:     %s\n', mat2str(size(expData)));
        return;
    end

    % Extract p-values safely
    p_Time = NaN; p_Group = NaN; p_Interaction = NaN;
    if numel(p) >= 4
        p_Time = p(1);
        p_Group = p(2);
        p_Interaction = p(4);
    elseif numel(p) == 2
        p_Time = p(1);
    end

    % Annotate top-left
    textX = 1;
    lineSpacing = 0.06 * yRange;
    yTop = yMax + 0.20 * yRange;  % Lower than the expanded ylim

    text(textX, yTop, sprintf('Group: p = %.3f', p_Group), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'Color', colorForP(p_Group), 'FontSize', 10);
    text(textX, yTop - lineSpacing, sprintf('Time: p = %.3f', p_Time), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'Color', colorForP(p_Time), 'FontSize', 10);
    text(textX, yTop - 2 * lineSpacing, sprintf('Interaction: p = %.3f', p_Interaction), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'Color', colorForP(p_Interaction), 'FontSize', 10);


    % Posthoc: Bonferroni-corrected t-tests if interaction is significant
    if ~isnan(p_Interaction) && p_Interaction < 0.05
        offset = 0.05 * range(allY);
        alpha = 0.05;
        m = numTimePoints;

        for t = 1:numTimePoints
            try
                [~, rawP] = ttest2(controlData(:, t), expData(:, t));
            catch
                rawP = NaN;
            end

            % Bonferroni correction
            pVal = min(rawP * m, 1);

            % Marker
            if pVal < 0.001
                marker = '***';
            elseif pVal < 0.01
                marker = '**';
            elseif pVal < 0.05
                marker = '*';
            else
                marker = '';
            end

            if ~isempty(marker)
                yPos = max(meanControl(t) + semControl(t), meanExp(t) + semExp(t)) + offset;
                text(xVals(t), yPos, marker, ...
                    'HorizontalAlignment', 'center', 'FontSize', 14, 'Color', colorForP(pVal));
            end
        end

        % Optional note
        text(textX, yTop - 3 * lineSpacing, '(Bonferroni-corrected posthoc)', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
            'Color', 'k', 'FontSize', 10);
    end
end

hold off;
end

% ======== Helper Function ========
function c = colorForP(p)
if isnan(p)
    c = 'k';
elseif p < 0.05
    c = 'r';
else
    c = 'k';
end
end
