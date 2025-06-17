function plotSegmentedData(Summary_ctrl, Summary_exp,  saveFigures, saveFormat, outputDir, enableStats)
% PLOTSEGMENTEDDATA Generates plots for segmented data and combined figures.
%   This function takes control and experimental data structures and the
%   sessions to include. It creates segmented plots and combined summary
%   figures for analysis. Optionally, figures can be saved to file.

% Default values if not provided
if nargin < 3 || isempty(saveFigures)
    saveFigures = 0;  % Default to not saving
end
if nargin < 4 || isempty(saveFormat)
    saveFormat = '.png';  % Default file format
end

if nargin < 5 || isempty(outputDir)
    outputDir = [];
end

if nargin < 6 || isempty(saveFormat)
    enableStats = 0;  
end


sessionNames = vertcat(Summary_ctrl.sessionNames) + " " + vertcat(Summary_ctrl.Sex);
%sessionNames = 'stim' + extractAfter(sessionNames,'stim');
sessionNames = string(regexp(sessionNames, 'day\d+_?(.*)', 'tokens', 'once'));
sessionIncluded = 1:size(Summary_ctrl,2);
region = Summary_ctrl.Region;
numSegments = size(Summary_ctrl(1).PortEntry_Seg,2);
fieldNames = fieldnames(Summary_ctrl);
segmentFieldNames = fieldNames(contains(fieldNames, '_Seg'));

%% Plot Averaged Segments
for fd = 1:length(segmentFieldNames)
    figure;
    %drawnow;
    set(gcf, 'Position', get(0, 'ScreenSize'));  % Make figure fullscreen
    im = 1;
    for session = sessionIncluded
        subplot(2, ceil(length(sessionIncluded)/2), im);
        [Nc, Ne] = plotTimeSeriesWithSEM_Multi(Summary_ctrl(session).(segmentFieldNames{fd}), Summary_exp(session).(segmentFieldNames{fd}), enableStats);
        if im == 1
            legend('show');
        end
        xlabel('session progression');
        xticks(1:numSegments);  % Tick positions

        % Generate percentage labels
        percentLabels = arrayfun(@(x) sprintf('%.0f%%', x), linspace(100/numSegments, 100, numSegments), 'UniformOutput', false);

        xticklabels(percentLabels);

        ylabel('Mean ± SEM');
        if ismissing(sessionNames(session))
            title(sprintf('sessionNames missing'));
        else
            title(sprintf('%s\nN = %d (C), %d (E)', sessionNames(session), Nc, Ne));
        end

        im = im + 1;
    end
    sgtitle("Segmented Plot: " + region + " " + extractBefore(segmentFieldNames{fd}, '_'));
    if saveFigures
        filename = region + "_Segmented_" + extractBefore(segmentFieldNames{fd}, '_') + saveFormat;
        filepath = fullfile(outputDir , filename);
        saveas(gcf, filepath);
    end
end

%% Plot Total Across Sessions

% Extract Field Names for Calculations
fieldNamesTotal = fieldNames(contains(fieldNames, '_Total'));
fieldNamesBouts = fieldNames(contains(fieldNames, '_SumBouts'));
fieldNamesFR = fieldNames(contains(fieldNames, '_SumFR'));

% Combine All Field Names
allFieldNames = [fieldNamesTotal; fieldNamesBouts; fieldNamesFR];
totalSubplots = length(allFieldNames);

% Create Figure and Determine Subplot Layout
figure;
set(gcf, 'Position', get(0, 'ScreenSize'));  % Make figure fullscreen
nRows = 3;
nCols = ceil(totalSubplots / nRows);

% Iterate Through All Field Names and Plot
subplotIdx = 1;
for fd = 1:length(allFieldNames)
    im = 1;
    for session = sessionIncluded
        AveC.(allFieldNames{fd}){im} = Summary_ctrl(session).(allFieldNames{fd});
        AveE.(allFieldNames{fd}){im} = Summary_exp(session).(allFieldNames{fd});
        im = im + 1;
    end

    subplot(nRows, nCols, subplotIdx);
    plotTimeSeriesWithSEM(AveC.(allFieldNames{fd}), AveE.(allFieldNames{fd}),enableStats);
    if fd == 1 
        legend('show');
    end
    xticklabels(sessionNames);
    ylabel('Mean ± SEM');
    title(extractBefore(allFieldNames{fd}, '_'));
    subplotIdx = subplotIdx + 1;
end
sgtitle('Total # Across Sessions: ' + region);
if saveFigures
    filename = region + "_Totals" + saveFormat;
    filepath = fullfile(outputDir , filename);
    saveas(gcf, filepath);
end
end
