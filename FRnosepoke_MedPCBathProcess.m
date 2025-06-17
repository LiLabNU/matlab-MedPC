clear all
close all
%% Step 1: Get folder directories and saved MedPC Mat filename
serverPath = "\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\PBS\LiPatel_Labs";
motherDir = "Personal_Folders\Valen\Optogenetics experiments combined\BNST";
motherDir = fullfile(serverPath, motherDir);
MedPCSaveName = 'processedMedPCData.mat';

datasetFolders = ["Females";...
    
    "Males",...
    ];

keywords = "FR3_"; % keywords for the folder you want to include

%% Step 2: Checks for previously processed MedPC data. If not, initiate the process.
% Full path to the file
processedDataFile = fullfile(motherDir, MedPCSaveName);
runProcess = false;

% Check if the .mat file exists
if exist(processedDataFile, 'file')
    % Load the file
    load(processedDataFile);
    fprintf('✅ Loaded preprocessed MedPC data from %s\n', processedDataFile);
    % Check if all dataset folders are included   
    loadedFolders = cellstr(datasetDir(:));
    expectedFolders = cellstr(fullfile(motherDir, datasetFolders(:)));    
    allIncluded = all(ismember(expectedFolders, loadedFolders));

    if allIncluded
        disp('✅ All expected dataset folders are present.');        
    else
        warning('⚠️ Some dataset folders are missing in the loaded datasetDir.');
        runProcess = true;
    end

else
    % File not found, set the flag
    warning('⚠️ Processed MedPC data not found in: %s', motherDir);
    runProcess = true;
end
% run the MedPC process loop
if  runProcess == true
    fprintf('Extracting MedPC data now...\n');

    datasetDir = fullfile(motherDir, datasetFolders);

    % Preallocate containers
    allData = cell(1, length(datasetDir));
    allTrialTS = cell(1, length(datasetDir));
    allTotals = cell(1, length(datasetDir));
    allDetails = cell(1, length(datasetDir));
    allAvg = cell(1, length(datasetDir));
    allMice = cell(1, length(datasetDir));
    allFolderName = cell(1, length(datasetDir));
    allCurrDir = cell(1, length(datasetDir));

    % Loop to process and store everything
    for f = 1:length(datasetDir)
        [data, trialTS, Totals, details, Avg, mice, folderName, currDir] = ...
            Nosepoke_MedPC_Valen_batchprocess(keywords, datasetDir(f));

        allData{f} = data;
        allTrialTS{f} = trialTS;
        allTotals{f} = Totals;
        allDetails{f} = details;
        allAvg{f} = Avg;
        allMice{f} = mice;
        allFolderName{f} = folderName;
        allCurrDir{f} = currDir;
    end

    % Save all the processed variables to as processedMedPCData.mat in the motherDir
    save(fullfile(motherDir, MedPCSaveName), ...
        'allData', 'allTrialTS', 'allTotals', 'allDetails', ...
        'allAvg', 'allMice', 'allFolderName', 'allCurrDir', 'datasetDir', '-v7.3');

    fprintf('Extracting MedPC data finished. Data saved as: %s\n', processedDataFile);
end

%% Step 3: Calculate Summary and Segments
refEvent = 'Sucrose'; % event name (in trialTS) that you want to segment your data based on
minSucEvent = 10; % exclude mice if they have less than X number of sucrose intakes


%
numSegments = 5; % how many segmentations do you want
minRefEvent = 10; % exclude mice if they have less than N number of selected refEvent

%
numSegTrials = 0;
minSegNum = 3;

SegmentSaveName = 'SummaryResults.mat';
MedPCBin = 100; % N number of bins for 1s

for f = 1:length(datasetDir)
    data = allData{f};
    trialTS = allTrialTS{f};
    Totals = allTotals{f};
    details = allDetails{f};
    Avg = allAvg{f};
    mice = allMice{f};
    folderName = allFolderName{f};
    currDir = allCurrDir{f};

    % get control and experimental mouse IDs
    temp = arrayfun(@(x) extractAfter(x, max(0, strlength(x) - 1)), mice, 'UniformOutput', false);
    temp = string(temp);
    nonEmptyColumns = find(all(temp ~= "", 1));

    exp =  mice(strcmp(temp(:,nonEmptyColumns(1)),'E'),1)';
    ctrl = mice(strcmp(temp(:,nonEmptyColumns(1)),'C'),1)';

    clear segmentedSum miceInclude
    for session = 1:size(data, 2)
        FRdetails_sess = cellfun(@(x) x, details.FRdetails(:,session), 'UniformOutput', false);
        PEdetails_sess = cellfun(@(x) x, details.PEdetails(:,session), 'UniformOutput', false);

        % segmenting into a fixed number of segmentations
        if numSegments ~=0 && numSegTrials ==0
            [segmentedSum(session,:), miceInclude(:,session)] = segmentedSumCalculation(data{session}.(refEvent), ...
                minSucEvent, minRefEvent, data{session}, MedPCBin, numSegments, trialTS{session}, ...
                FRdetails_sess, PEdetails_sess, refEvent, mice(:,session));

        % segmenting based on a fixed number of trials
        elseif numSegments ==0 && numSegTrials ~=0    
            [segmentedSum(session,:), miceInclude(:,session)] = segmentEdges(data{session}.(refEvent), minSucEvent, minSegNum, data{session}, trialTS{session}, ...
                FRdetails_sess, PEdetails_sess, numSegTrials, mice);

        else
            warning('Please check numSegments and numSegTrials. One of them has to be zero in order to specify a segmentation method!');
        end

    end

    % combine summaries
    if f == 1
        [Summary_ctrl, Summary_exp] = summarizeSessions(miceInclude, ctrl, exp, Totals, Avg, segmentedSum, folderName, currDir);
    else
        [newSummary_ctrl, newSummary_exp] = summarizeSessions(miceInclude, ctrl, exp, Totals, Avg, segmentedSum, folderName, currDir);

        Summary_ctrl = [Summary_ctrl, newSummary_ctrl];
        Summary_exp = [Summary_exp, newSummary_exp];
    end
end

% save data
save(fullfile(motherDir, SegmentSaveName), 'Summary_ctrl', 'Summary_exp');
disp("MedPC data Summary and Segments done.");
%% Step 4: plot selected sessions from selected dataset and save the figures in the mother folder
saveFigures = 1;
saveFormat = '.png';
enableStats = 0;

% Nucleus Accumbens
region = 'April_Opto_inhibition_CeA';
sessionIdx = contains(vertcat(Summary_ctrl.sessionNames), "stim") | contains(vertcat(Summary_ctrl.sessionNames), "inhib");
regionIdx = strcmpi(vertcat(Summary_ctrl.Region), region);
plotIdx = regionIdx & sessionIdx;
if sum(plotIdx) ~=0
    plotSegmentedData(Summary_ctrl(plotIdx), Summary_exp(plotIdx), saveFigures, saveFormat, motherDir, enableStats);
end

% Central Amygdala
region = 'Central Amygdala';
sessionIdx = contains(vertcat(Summary_ctrl.sessionNames), "stim");
regionIdx = strcmpi(vertcat(Summary_ctrl.Region), region);
plotIdx = regionIdx & sessionIdx;
if sum(plotIdx) ~=0
    plotSegmentedData(Summary_ctrl(plotIdx), Summary_exp(plotIdx), saveFigures, saveFormat, motherDir, enableStats);
end

% BNST
region = 'BNST';
sessionIdx = contains(vertcat(Summary_ctrl.sessionNames), "stim");
regionIdx = strcmpi(vertcat(Summary_ctrl.Region), region);
plotIdx = regionIdx & sessionIdx;
if sum(plotIdx) ~=0
    plotSegmentedData(Summary_ctrl(plotIdx), Summary_exp(plotIdx), saveFigures, saveFormat, motherDir, enableStats);
end

%close all




















