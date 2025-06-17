  clear all
%close all
%% get files and port entries
% motherFolder = 'R:\PBS\LiPatel_Labs\Personal_Folders\Ziqing\SA\SA-pilot3\MedPC_MedPC_data';
% [Folders, fileNames] = getDirContents(motherFolder, includeSubfolders);

[MedPCfile, folderDir] = GetFilesFromFolder(0,'.txt');
[MedPC_data, MedPC_trialTS, AnimalIDcell] = SucroseSA_MedPC2mat(MedPCfile, folderDir);

%% sort results based on animal's ID
mice = AnimalIDcell(:,1);
[miceID, idx] = sort(mice);
AnimalIDcell = AnimalIDcell(idx,:);

fields = fieldnames(MedPC_data);
for fd = 1:length(fields)
    field = fields{fd};
    MedPC_data.(field)= MedPC_data.(field)(idx,:);
end

fields = fieldnames(MedPC_trialTS);
for fd = 1:length(fields)
    field = fields{fd};
    MedPC_trialTS.(field)= MedPC_trialTS.(field)(:,idx);
end

%% raster plot for individual mice
bin = 60; % plotting 60s of MedPC_data per line in the subplot
imColumns = 6; % define numbers of columns of the subplot
imRows = ceil(numel(MedPCfile)/imColumns); % define numbers of rows of the subplot
legendLabels = {'Reward', 'ActiveNP', 'InactiveNP', 'RewardTrigNP', 'PortOccupancy'};

Plot_individualRaster_SA(MedPC_data, AnimalIDcell, [imRows imColumns], bin, legendLabels);

%% Plot port entry PSTHs for sucrose
pre = 10;
post = 20;
MedPCBin = 100; % 100 bins = 1s if bin size is 10ms

event = 'RewardTrigNP';
ts = MedPC_trialTS.(event);
psth_matrix = {};
for i = 1:numel(ts)
    onsets = ts{i};
    for j = 1:numel(onsets)
        startIdx = onsets(j) - pre * MedPCBin;
        endIdx = onsets(j) + post * MedPCBin;

        % Determine the expected trial length
        expectedLength = pre * MedPCBin + post * MedPCBin + 1;

        % Handle negative start index
        if startIdx < 1
            padStart = abs(startIdx) + 1; % Calculate how many zeros are needed
            startIdx = 1; % Adjust start index to valid range
        else
            padStart = 0; % No need to pad if start index is valid
        end

        % Handle end index exceeding available MedPC_data
        if endIdx > size(MedPC_data.PortEntry, 2)
            padEnd = endIdx - size(MedPC_data.PortEntry, 2);
            endIdx = size(MedPC_data.PortEntry, 2); % Adjust end index to valid range
        else
            padEnd = 0; % No need to pad if end index is valid
        end

        % Extract the valid portion of the MedPC_data
        validData = MedPC_data.PortEntry(i, startIdx:endIdx);

        % Pad with zeros if needed
        paddedData = [zeros(1, padStart), validData, zeros(1, padEnd)];

        % Ensure consistent trial length
        if length(paddedData) < expectedLength
            paddedData = [paddedData, zeros(1, expectedLength - length(paddedData))];
        end

        psth_matrix{i}(j, :) = paddedData;
    end
end

figure
imColumns = 6; % define numbers of columns of the subplot
imRows = ceil(numel(MedPCfile)/imColumns); % define numbers of rows of the subplot
binFactor = 10;
tickSpace = 4;
for i = 1:numel(psth_matrix)
    subplot(imRows,imColumns,i)
    binwindow = size(psth_matrix{i},2)/binFactor;
    plotAverageWithSEM(imresize(psth_matrix{i},[size(psth_matrix{i},1),binwindow],'bilinear'));

    % Dynamically set x ticks evenly spaced
    xticks_positions = linspace(1, binwindow, tickSpace);  % 5 evenly spaced ticks
    xticks(xticks_positions);
    % Calculate corresponding time labels
    time_labels = linspace(-pre, post, tickSpace);
    xticklabels(round(time_labels, 1)); % Round to avoid long decimals
    xlabel('Time from Tone onset (s)')
    ylabel('Port Entry Probability')
    title(miceID(i))
end


%% Calculate # of events
for i=1:length(MedPCfile) % for i equal to 1 to the # of medPC files (i.e., animals)
    Totals.PortOccupancy(i,1) = sum(MedPC_data.PortEntry(i,:))./MedPCBin; % port occupancy in seconds
    Totals.Reward(i,1) = sum(MedPC_data.Reward(i,:));
    Totals.ActiveNP(i,1) = sum(MedPC_data.ActiveNP(i,:));
    Totals.InactiveNP(i,1) = sum(MedPC_data.InactiveNP(i,:));
    Totals.RewardTrigNP(i,1) = sum(MedPC_data.RewardTrigNP(i,:));
end

%%
psth_matrix = psth_matrix';
clearvars -except MedPC_data Totals MedPC_trialTS miceID psth_matrix


%%% code below is not tested
% %% Calculate port entry related readouts
% MedPCBin = 100; % 100 bin is 1s for medpc if binsize is 10ms
% fr = 1;
% clear PEdetials SumPE FRdetails SumFR
% for i=1:length(MedPCfile)
%     eventTimestamps = trialTS.Sucrose{i};
%     if ~isempty(eventTimestamps)
%         [PEdetials{i}, SumPE(i)] = findBouts(MedPC_data.PortEntry(i,:), MedPCBin, eventTimestamps, AnimalIDcell(i,1), 'mean');
%         [FRdetails{i}, SumFR(i)] = findFRduration(trialTS.FirstActive{i}, eventTimestamps, trialTS.ActiveNP{i},fr);
%     end
% end
% 
% %% Segment timeseries
% clear segmentedSum
% numSegments = 5; % how many segmentations do you want
% minEvent = 10;
% refEvent = 'Sucrose';
% refData = MedPC_data.(refEvent); % do you want to segment based on Sucrose deliver?
% %TimeSeries = ones([1 length(MedPC_data.Sucrose)]);
% 
% % find position where first 25% 50% 75% last sucrose delivery is completed
% Pos = ceil(findPositions(refData,minEvent)/MedPCBin);
% 
% fieldNames = fieldnames(MedPC_data);
% numFields = length(fieldNames);
% % Loop through each file
% for i = 1:length(MedPCfile)
%     RefSeries = refData(i,:);
%     if sum(RefSeries) >= minEvent
%         tempResults = calculateSegment(MedPC_data, RefSeries, numSegments, i, 'event');  % Assume this function is adjusted
%         for k = 1:numFields
%             fieldName = fieldNames{k};
%             segmentedSum.(fieldName)(i, :) = tempResults.(fieldName);
%         end
%     else
%         for k = 1:numFields
%             fieldName = fieldNames{k};
%             segmentedSum.(fieldName)(i, :) = 0;
%         end
%     end
% end
% 
% % Assign matrices to base workspace
% for k = 1:numFields
%     assignin('base', fieldNames{k}, segmentedSum.(fieldNames{k}));
% end
% 
% clear averages
% % Loop through each file
% for i = 1:length(MedPCfile)
%     RefSeries = refData(i,:);
%     TimeSeries = trialTS.(refEvent){i};
%     if ~isempty(FRdetails{i})
%         MedPC_dataSeries = FRdetails{i}.FRDuration + PEdetials{i}.firstBoutLatency;
%         if sum(RefSeries) >= minEvent
%             averages(i,:) = segmentTimestamps(TimeSeries, PEdetials{i}.eventTime, MedPC_dataSeries, numSegments,'event');
%         else
%             averages(i,:) = 0;
%         end
%     end
% end
% averages = ceil(averages/MedPCBin);
% segmentedSum.('LatencyNPtoPE') = averages;
% 
% 
% clear averages
% % Loop through each file
% for i = 1:length(MedPCfile)
%     RefSeries = refData(i,:);
%     TimeSeries = trialTS.(refEvent){i};
%     if ~isempty(FRdetails{i})
%         MedPC_dataSeries = FRdetails{i}.FRcompletion;
%         if sum(RefSeries) >= minEvent
%             averages(i,:) = segmentTimestamps(TimeSeries, trialTS.FirstActive{i}, MedPC_dataSeries, numSegments,'event');
%         else
%             averages(i,:) = 0;
%         end
%     end
% end
% segmentedSum.('FRCompletion') = averages;
% 
% clear averages
% % Loop through each file
% for i = 1:length(MedPCfile)
%     RefSeries = refData(i,:);
%     TimeSeries = trialTS.(refEvent){i};
%     if ~isempty(FRdetails{i})
%         MedPC_dataSeries = PEdetials{i}.interEventITI;
%         if sum(RefSeries) >= minEvent
%             averages(i,:) = segmentTimestamps(TimeSeries, PEdetials{i}.eventTime, MedPC_dataSeries, numSegments,'event');
%         else
%             averages(i,:) = 0;
%         end
%     end
% end
% %averages(averages==0) = 360000/numSegments;
% averages = ceil(averages/MedPCBin);
% segmentedSum.('interEventITI') = averages;
% 
% group1 = ["0";"11263";"11265";"7974", "7980"]; % Males control
% group2 = ["1";"11264";"7973";"7975", "7981", "7982"]; % Females CeA opto; ChR2 bad
% 
% 
% temp = segmentedSum.ActiveNP;
% temp = segmentedSum.LatencyNPtoPE;
% temp = segmentedSum.FRCompletion;
% miceToUse = ismember(mice,[group1;group4]);
% tempE = temp(miceToUse,:);
% miceToUse = ismember(mice,[group3;group6]);
% tempC = temp(miceToUse,:);
% 
% miceToUse = ismember(mice,[group7;group9]);
% tempE = temp(miceToUse,:);
% miceToUse = ismember(mice,[group8;group10]);
% tempC = temp(miceToUse,:);
% 
% miceToUse = ismember(mice,[group13]);
% tempE = temp(miceToUse,:);
% miceToUse = ismember(mice,[group12]);
% tempC = temp(miceToUse,:);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% for i = 1: size(results,2)
%     if ~isempty(results{i})
%         firstBoutLatency{i} = results{i}.firstBoutLatency;
%         temp = find(nanzscore(firstBoutLatency{i})>2,1,'first');
%         %temp = find(firstBoutLatency{i}(6:end)>(nanmean(firstBoutLatency{i})+nanstd(firstBoutLatency{i})*2),1,"first");
%         if ~isempty(temp)
%             idx(i) = temp+5;
%         else
%             idx(i) = length(firstBoutLatency{i});
%         end
%     else
%         firstBoutLatency{i} = [];
%         idx(i) = 0;
%     end
% end
% 
% shk = ceil(idx/5);
% 
% control = [1:4 8:12]; % control idx for April Opto CeA Male
% %control = [4:10]; % control idx for April Opto CeA female
% control = [1 2 5 6 9 10 13 14]; % control idx for CRISPR PVT-CeA female
% 
% for i = 1: size(results,2)
%     if ~isempty(results{i}) & Totals.Sucrose_n(i) > 10
%         firstBoutLatency{i} = results{i}.firstBoutLatency;
%         if ismember(i,control) % [1:4 12:14] are the idx for control mice
%             color = [0, 0, 1, 0.5]; % control - blue
%         else
%             color = [1, 0, 0, 0.5]; % experimental - red
%         end
%         temp = firstBoutLatency{i};
%         temp = temp(~isnan(temp));
%         MedPC_dataCumu = generateCumulativePlots(temp','probability', color);
%         hold on
% 
%         idx(i) = find((MedPC_dataCumu)>0.5,1,'first');
%     end
% end
% 
% 
% 
% 
% %% animal included
% include = []; % for CRISPR PVT-CeA female
% 
% include = [9:14]; % for April Opto CeA female CONTROL
% include = [1:4]; % for April Opto CeA female EXPERIMENTAL
% 
% include = [1:4,12:14]; % for April Opto CeA male CONTROL
% include = [6:8]; % for April Opto CeA male EXPERIMENTAL
% 
% %% separating groups
% 
% % control group index
% control = [1 2 5 6 9 10 13 14]; % control idx for CRISPR PVT-CeA female
% control = [1:4 12:14]; % control idx for April Opto CeA Male
% control = [9:15]; % control idx for April Opto CeA female
% MedPC_dataSum = struct();
% a = 1;
% b = 1;
% fieldName = 'meanFirstBoutLatency';
% for i = 1: size(Sum,2)
%     if ismember(i,control)
%         if ~isempty(Sum(i).(fieldName))
%             MedPC_dataSum.ctrl(a) = Sum(i).(fieldName);
%         else
%             MedPC_dataSum.ctrl(a) = 0;
%         end
%         a = a+1;
%     else
%         if ~isempty(Sum(i).(fieldName))
%             MedPC_dataSum.exp(b) = Sum(i).(fieldName);
%         else
%             MedPC_dataSum.exp(b) = 0;
%         end
%         b = b+1;
%     end
% end
% 
% %%
% 
% % shock intensity schedule
% startIntensity = 0.02;
% multiply = 1.5;
% step = 5;
% % calculate shock intensity per sucrose delivery
% for i=1:length(MedPCfile)
%     temp(i) = sum(FRdetails{i}.FRcompletion);
% end
% maxV = max(temp);
% maxV = maxV(1);
% idx = find(temp==maxV);
% idx = idx(1);
% 
% %% plot shock cumulative
% figure
% for i=1:length(MedPCfile)
%     shkvec{i} = createShockVector(startIntensity, length(trialTS.Sucrose{i}), multiply, step);
%     if ismember(i,control) % [1:4 12:14] are the idx for control mice
%         color = [0, 0, 1, 0.5]; % control - blue
%     else
%         color = [1, 0, 0, 0.5]; % experimental - red
%     end
%     cumulatives{i} = generateCumulativePlots(FRdetails{i}.FRcompletion','raw',color);
%     hold on
% 
%     % simple linear regression
%     lm = fitlm(1:length(cumulatives{i}),cumulatives{i});
%     slope(i) = lm.Coefficients.Estimate(2);
% 
%     % find EC50
%     EC50_shock(i) = shkvec{i}(floor(length(shkvec{i})/2));
%     nm = 1;
%     for j = 1:length(FRdetails{i}.FRcompletion)
%         temp = FRdetails{i}.FRcompletion(j);
%         if temp == 1
%             shk{i}(j) = shkvec{i}(nm);
%             nm = nm+1;
%         else
%             if nm == 1
%                 shk{i}(j) = shkvec{i}(nm);
%             else
%                 shk{i}(j) = shkvec{i}(nm-1);
%             end
%         end
%     end
%     EC50_pos(i,:) = find(shk{i}==EC50_shock(i),1,'first');
% end
% 
% %% plot unfinished attempts cumulative
% figure
% for i=1:length(MedPCfile)
%     temp = FRdetails{i}.FRcompletion;
%     if temp(end) == 0
%         temp(end+1) = 1;
%     end
%     a = find(temp == 1);
%     numZeros = diff(a)-1;
%     unfinishedAtte{i} = numZeros;
%     if ismember(i,control) % the idx for control mice
%         color = [0, 0, 1, 0.5]; % control - blue
%     else
%         color = [1, 0, 0, 0.5]; % experimental - red
%     end
%     generateCumulativePlots(numZeros','raw',color);
%     hold on
% end
% ylabel('Unfinished Attempts')
% xticks(1:5:maxV)
% xticklabels(round(shkvec{idx}(1:5:end),2))
% xlabel('Shock Intensity (mA)')
% 
% %% calculate breakpoint
% timeout = 2; % if 2 consecutive 60s timeout, then we think it is the breakpoint.
% for i=1:length(MedPCfile)
%     bp = find(unfinishedAtte{i}>=timeout);
%     % ignore the first 10 trials
%     bp = bp(find(bp>10, 1,'first'));
%     if isempty(bp)
%         breakpoint(i,:) = max(shkvec{i}(1:end-1));
%     else
%         breakpoint(i,:) = shkvec{i}(bp-1);
%     end
%     maxShock(i,:) = max(shkvec{i}(1:end-1));
% end
% 
% %%
% for i=1:length(MedPCfile)
%     x = shkvec{i}';
%     y = results{i}.interEventBoutDuration;
%     R = corrcoef(x, y);
%     if length(R) == 1
%         coeff(i) = 0;
%     else
%         coeff(i) = R(1,2);
%     end
% end
% 




% %% calculate and plot trial average for cued
% win_time = 10;
% win_half = win_time*100; % 30s before and after cue onset
% binsize = 0.5; % 0.5s per bin
% TS = trialTS.Cue;
% events = ["ActiveNP"];
% for j = 1:length(events)
%     event = cell2mat(events(j));
%     MedPC_dataToCalculate = MedPC_data.(event);
%     if TS{1} ~= 0
%         trial_average.(event) = processPE(TS, win_half, binsize, MedPC_dataToCalculate);
%
%         % plot individual trial average
%         session = 1;
%         mice = 14;
%         plotcoloum = ceil(mice/1);
%
%         MedPC_dataToplot = trial_average.(event).mxt.rew_ms;
%         counter = 1;
%         for i = ((session-1)*mice)+1:(session)*mice
%             subplot(plotcoloum, ceil(mice/plotcoloum), counter)
%             plot(MedPC_dataToplot(i,:,:))
%             hold on
%             title(AnimalIDcell(i,1))
%             %ylim([-5 20])
%             ylabel('Probability')
%             xlim([0 size(MedPC_dataToplot,2)])
%             xticks([0 floor(size(MedPC_dataToplot,2)/2) size(MedPC_dataToplot,2)])
%             xticklabels([win_time/-1 0 win_time])
%             xlabel('Time from CS onset')
%             if counter == 1
%                 legend(event,'Location','northwest')
%             end
%             counter = counter +1;
%         end
%     end
% end
%
% %% some more complicated analysis
%
% % basic calculations
%  [beforeCueAnalysis, afterCueAnalysis] = analyzeNosepoke(trialTS);
% % average within each trial
% for i = 1:size(afterCueAnalysis,2)
%     meanStruct{i} = calculateMeansForStructArray(afterCueAnalysis{i});
% end
% % average every 5 trials that have the same shock intensity
% for i = 1:size(afterCueAnalysis,2)
%     reducedStruct{i} = averageFieldsInStruct(meanStruct{i}, 5);
% end
%
% % Define the path to the Excel file that has box information
% filename = 'Z:\PBS\LiPatel_Labs\Personal_Folders\Valen\01162024_NT_CRISPR_Females_VO\FR3_cue_shock_day18\Test1\MouseID_boxes.xlsx';
% % Read the entire Excel file
% MedPC_dataTable = readtable(filename);
% % Reorder the box information to match the order in AnimalIDcell
% tableIDcell = string(table2cell(MedPC_dataTable(:,1)));
% [~, loc] = ismember(tableIDcell,AnimalIDcell(:,1));
% MedPC_dataTable = MedPC_dataTable(loc,:);
% box = string(table2cell(MedPC_dataTable(:,2)));
% condition = string(table2cell(MedPC_dataTable(:,3)));
%
% % exclude any mice
% exclude = strcmp(box,"1") | strcmp(box,"3");
% condition = condition(~exclude,:);
% mouseID = AnimalIDcell(~exclude,1);
% %MedPC_dataToUse = meanStruct;
% MedPC_dataToUse = reducedStruct;
% MedPC_dataToUse = MedPC_dataToUse(~exclude);
%
% % Determine a suitable subplot grid size (for simplicity, aiming for a square grid)
% fieldNames = fieldnames(MedPC_dataToUse{1});
% gridSize = ceil(sqrt(length(fieldNames)));
% figure; % Create a new figure for the subplots
% result = struct();
% for i = 1:length(fieldNames)
%     MedPC_data2 = [];
%     fieldName = fieldNames{i};
%     for m = 1: size(MedPC_dataToUse,2)
%         current = MedPC_dataToUse{m};
%         MedPC_data = [];
%         for j = 1:length(current)
%             temp = current(j).(fieldName);
%             MedPC_data = [MedPC_data temp];
%         end
%         MedPC_data2(:,m) = MedPC_data;
%     end
%     result.(fieldName) = MedPC_data2;
%
%     %plot
%     subplot(gridSize, ceil(length(fieldNames)/gridSize), i);
%     errorbar_pn_hao(MedPC_data2(:,condition=="control")', [0, 0.4470, 0.7410]) % blue
%     hold on
%     errorbar_pn_hao(MedPC_data2(:,condition=="NT")', [0.8500, 0.3250, 0.0980]) % orange
%     hold off
%     title(sprintf(fieldName))
% end
%
% %%
% MedPC_dataToUse = trial_average.mxt.rew_ms;
% groups = [1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2];
%
% % Separate MedPC_data into two groups
% group1Data = MedPC_dataToUse(groups == 1, :);
% group2Data = MedPC_dataToUse(groups == 2, :);
%
% % Calculate averages for bins 1 to 60 (before) and bins 61 to 121 (after) for Group 1
% averageBefore60_Group1 = mean(group1Data(:, 1:60), 2);
% averageAfter60_Group1 = mean(group1Data(:, 61:end), 2);
%
% % Combine the averages for easier interpretation for Group 1
% averages_Group1 = [averageBefore60_Group1, averageAfter60_Group1];
%
% % Calculate averages for bins 1 to 60 (before) and bins 61 to 121 (after) for Group 2
% averageBefore60_Group2 = mean(group2Data(:, 1:60), 2);
% averageAfter60_Group2 = mean(group2Data(:, 61:end), 2);
%
% % Combine the averages for easier interpretation for Group 2
% averages_Group2 = [averageBefore60_Group2, averageAfter60_Group2];
%
%
% % Plot the results
% HaoBarErrorbar(averageBefore60_Group1, averageAfter60_Group1, [], 'mean');
% set(gca, 'XTickLabel', {'PreCue', 'PostCue'});
% title('Control')
% HaoBarErrorbar(averageBefore60_Group2, averageAfter60_Group2, [], 'mean');
% title('CRISPR')
% set(gca, 'XTickLabel', {'PreCue', 'PostCue'});
%
% HaoBarErrorbar(averageBefore60_Group1, averageBefore60_Group2, [], 'mean');
% set(gca, 'XTickLabel', {'Control', 'CRISPR'});
% title('PreCue')
% HaoBarErrorbar(averageAfter60_Group1, averageAfter60_Group2, [], 'mean');
% title('PostCue')
% set(gca, 'XTickLabel', {'Control', 'CRISPR'});
%
% HaoBarErrorbar(averageAfter60_Group1-averageBefore60_Group1, averageAfter60_Group2-averageBefore60_Group2, [], 'mean');
% set(gca, 'XTickLabel', {'Control', 'CRISPR'});
% title('PostCue - PreCue')
%
% %% converting trialTS to match SLEAP prediction
%
% % define frame rate of behavioral videos
% framerate = 15;
%
% % for day 18
% indices = [5 13	6	14	7	15	8	16	1	9	2	10	3	11	4	12];
% mouseID = AnimalIDcell(indices,1);
%
% % Make a copy of the original structure to the new structure 'trialTS_slp'
% trialTS_slp = trialTS;
%
% % Get a list of all field names in the original structure
% fieldNames = fieldnames(trialTS);
%
% % Loop through each field name in the original structure
% for i = 1:numel(fieldNames)
%     % Current field name
%     currentFieldName = fieldNames{i};
%
%     % Check if the field contains a cell array in the original structure
%     if iscell(trialTS.(currentFieldName))
%         % Get the current cell array from the original structure
%         currentCellArray = trialTS.(currentFieldName);
%
%         % Initialize a new cell array for the rearranged and modified MedPC_data
%         newCellArray = cell(size(currentCellArray));
%
%         % Loop through each cell in the current field's cell array to apply modifications
%         for j = 1:numel(currentCellArray)
%             % Apply the operation (divide by 100 and then multiply by 15) to each cell's value
%             % Ensure the cell contains numeric MedPC_data
%             if isnumeric(currentCellArray{j})
%                 newCellArray{j} = currentCellArray{j} / 100 * 15;
%             else
%                 warning('Non-numeric MedPC_data encountered in field %s, cell %d. Skipping...', currentFieldName, j);
%             end
%         end
%
%         % Rearrange the cells in the new cell array according to the provided indices
%         % Ensure that the length of indices matches the number of cells
%         if length(indices) == numel(currentCellArray)
%             rearrangedCellArray = cell(size(newCellArray));
%             for k = 1:numel(newCellArray)
%                 rearrangedCellArray{k} = newCellArray{indices(k)};
%             end
%             % Update the field in the new structure with the rearranged cell array
%             trialTS_slp.(currentFieldName) = rearrangedCellArray;
%         else
%             warning('Mismatch in the number of indices and cells in field %s. Skipping rearrangement...', currentFieldName);
%         end
%     else
%         warning('Field %s does not contain a cell array. Skipping...', currentFieldName);
%     end
% end
%

