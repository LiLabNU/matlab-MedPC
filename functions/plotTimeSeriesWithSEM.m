function plotTimeSeriesWithSEM(controlData, expData, enableStats)
    % Function to plot control and experimental data across time points with SEM error bars
    % Inputs:
    %   - controlData: Cell array where each cell contains a vector of control samples for a time point
    %   - expData: Cell array where each cell contains a vector of experimental samples for a time point
    %   - enableStats: Logical flag; if true, perform statistical tests (unpaired t-test) at each time point
    %
    % Note: This function assumes that each cell in controlData and expData contains numeric data.
    
    % Number of time points
    numTimePoints = length(controlData);

    % Preallocate for means and SEMs
    meanControl = zeros(1, numTimePoints);
    semControl = zeros(1, numTimePoints);
    meanExp = zeros(1, numTimePoints);
    semExp = zeros(1, numTimePoints);

    % Calculate mean and SEM for each time point
    for t = 1:numTimePoints
        % Control group
        meanControl(t) = mean(controlData{t}, 'omitnan');
        semControl(t) = std(controlData{t}, 'omitnan') / sqrt(length(controlData{t}));

        % Experimental group
        meanExp(t) = mean(expData{t}, 'omitnan');
        semExp(t) = std(expData{t}, 'omitnan') / sqrt(length(expData{t}));
    end

    % Define the time points (x-axis)
    timePoints = 1:numTimePoints;

    % Plot the data
    
    hold on;

    % Plot control data with error bars
    errorbar(timePoints, meanControl, semControl, '-o', 'LineWidth', 1.5, 'Color', 'b', 'DisplayName', 'Ctrl');

    % Plot experimental data with error bars
    errorbar(timePoints, meanExp, semExp, '-o', 'LineWidth', 1.5, 'Color', 'r', 'DisplayName', 'Exp');

    % Customize plot appearance
    xticks(timePoints);
    grid on;
    
    if enableStats
        % Determine an offset for significance markers (5% of the overall range)
        yAll = [meanControl + semControl, meanExp + semExp];
        offset = 0.05 * (max(yAll) - min(yAll));

        for t = 1:numTimePoints
            % Use curly braces to extract numeric vectors from the cell arrays.
            [~, p] = ttest2(controlData{t}, expData{t});

            % Determine significance marker based on p-value thresholds
            if p < 0.001
                sigMarker = '***';
            elseif p < 0.01
                sigMarker = '**';
            elseif p < 0.05
                sigMarker = '*';
            else
                sigMarker = '';
            end

            % If significant, place the marker above the highest error bar at this time point.
            if ~isempty(sigMarker)
                yPos = max(meanControl(t) + semControl(t), meanExp(t) + semExp(t)) + offset;
                text(timePoints(t), yPos, sigMarker, 'HorizontalAlignment', 'center', ...
                    'FontSize', 14, 'Color', 'k');
            end
        end
    end

    hold off;
end
