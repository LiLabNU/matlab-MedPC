% add PBS\LiPatel_Labs\General\MatlabCode folders to path
% add my DLC folder to path (for stderr code)

clear all
close all

%% get files and port entries
[MedPCfile, folderDir] = GetFilesFromFolder(0,'.txt');
MedPCfile = MedPCfile(~startsWith(MedPCfile, '~$')); % delete temp/lock files
% MedPCfile = MedPCfile(contains(MedPCfile,'!'));
% folderDir = folderDir(contains(MedPCfile,'!'));

[allport, ~ , AnimalIDcell,trialtype,trialTS] = PortEntry_MedPC2mat(MedPCfile, folderDir);

%% get psth
win_time = 15; % how many seconds 
win_half = win_time*100; % before and after cue onset
cue_dur = 5; % how long is the cue !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! edit
binsize = 0.5; % 0.5s per bin
trialtypes = unique(trialtype{1}); % get trial types
for i = 1:length(trialtypes)
    ty = trialtypes(i);
    cue_onset = [trialtype{i}, trialTS{i}];
    pe(i) = processPE(cue_onset, win_half, cue_dur, binsize, allport, ty);
end
%save('port entry.mat','pe','trialNames', 'AnimalIDcell','allport');


%% plot individual 
session = 1; % how many sessions you have
mice = size(AnimalIDcell,1); % how many mice you have

trialNames = ["Sucrose CS"; "Shock CS"; "Neutral CS"];
trialToPlot = [1 2 4]; % which trial type(s) you want to plot: 1 is sucrose, 2 is shock, 4 is neutral !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! edit

figure;
for j = 1: length(trialToPlot)
    ty = find(trialToPlot(j)==trialtypes);
    dataToplot = pe(ty).mxt.rew_ms; % this is plotting the trial average of the selected trial type
    counter = 1;
    for i = ((session-1)*mice)+1:(session)*mice
        plotcoloum = ceil(mice/3);
        subplot(plotcoloum, ceil(mice/plotcoloum), counter)
        plot(dataToplot(i,:))
        hold on
        title(AnimalIDcell(i,1))
        ylim([0 1])
        ylabel('P(port entry)')
        xlim([0 size(dataToplot,2)])
        xticks([0 floor(size(dataToplot,2)/2) size(dataToplot,2)])
        xticklabels([win_time/-1 0 win_time])
        xlabel('Time from CS onset (s)')
        warning off
        legend(trialNames,'Location','northwest')
        warning on
        counter = counter +1;
    end
end

day = regexp(folderDir{1}, 'Day\d+', 'match', 'once');
fileName = "Plot_" + day + ".fig";
figPath = fullfile(folderDir{1}, fileName);
savefig(gcf, figPath);

%{
figure
for j = 1: length(trialToPlot)
    ty = find(trialToPlot(j)==trialtypes);
    dataToplot = pe(ty).mxt.rew_ms; % this is plotting the trial average of the selected trial type
    subplot(1, 3, j)
    errorbar_pn_hao(dataToplot);
    title(trialNames(j))
    ylim([0 1])
    ylabel('P(port entry)')
    xlim([0 size(dataToplot,2)])
    xticks([0 floor(size(dataToplot,2)/2) size(dataToplot,2)])
    xticklabels([win_time/-1 0 win_time])
    xlabel('Time from CS onset (s)')
end
%}

%% get blocks structures (each individual)
blkName = ["Appetitive1";"Aversive1";"Appetitive2";"Aversive2"];
blk = {1:30; 31:60; 61:90; 91:120};
for i = 1:length(trialtypes)
    ty = trialtypes(i);
    for b = 1:numel(blkName)
        cue_onset = [trialtype{i}(blk{b}), trialTS{i}(blk{b})];
        peBLK(i).(blkName(b)) = processPE(cue_onset, win_half, cue_dur, binsize, allport, ty);
    end
end

%% get blocks structures (average blocks)
blkName = ["AppetitiveAve";"AversiveAve"];
blk = {[1:30 61:90]; [31:60 91:120]};
for i = 1:length(trialtypes)
    ty = trialtypes(i);
    for b = 1:numel(blkName)
        cue_onset = [trialtype{i}(blk{b}), trialTS{i}(blk{b})];
        peBLK(i).(blkName(b)) = processPE(cue_onset, win_half, cue_dur, binsize, allport, ty);
    end
end

%% Plot individual for blocks (peBLK)
fields = fieldnames(peBLK(1));  % get block field names from the first trial type

for f = 1:length(fields)
    blockField = fields{f};       % e.g. 'Appetitive1', 'Aversive1', etc.
    figure; % Create a new figure for this block field
    for j = 1:length(trialToPlot)
        trialIdx = trialToPlot(j);
        ty = find(trialIdx == trialtypes); % Find the index in trialtypes corresponding to the desired trial type value
        if isempty(ty)
            continue; % skip if this trial type is not present
        end
        dataToplot = peBLK(ty).(blockField).mxt.rew_ms; % Get the data for this block field and trial type
        counter = 1;
        for i = ((session-1)*mice)+1:(session*mice)
            plotCols = ceil(mice/3); % Determine the number of subplot columns (adjust as needed)
            subplot(plotCols, ceil(mice/plotCols), counter)
            plot(dataToplot(i,:))
            hold on
            title(AnimalIDcell(i,1))
            ylim([0 1])
            ylabel('P(port entry)')
            xlim([0 size(dataToplot,2)])
            xticks([0 floor(size(dataToplot,2)/2) size(dataToplot,2)])
            xticklabels([win_time/-1 0 win_time])
            xlabel('Time from CS onset (s)')
            warning off
            legend(trialNames, 'Location', 'northwest')
            warning on
            counter = counter + 1;
        end
        % Add an overall title for the figure showing block and trial type info
        sgtitle(sprintf('%s', blockField))
    end
    figName = "PE_" + blockField + ".fig";
    figPath = fullfile(folderDir{1}, figName);
    savefig(gcf, figPath);
end

%% Group average plotting for blocks (peBLK)
fields = fieldnames(peBLK(1));  % Get block field names (assumed to be the same for each trial type)

for f = 1:length(fields)
    blockField = fields{f};  % e.g., 'Appetitive1', 'Aversive1', etc.
    
    % Create a new figure for the current block field
    figure;
    for j = 1:length(trialToPlot)
        % Identify the trial type index corresponding to the current trialToPlot value
        ty = find(trialToPlot(j) == trialtypes);
        if isempty(ty)
            continue;  % Skip if the trial type is not found
        end
        
        % Extract the data for group average plotting from the block structure
        dataToplot = peBLK(ty).(blockField).mxt.rew_ms;
        
        % Create a subplot for the current trial type
        subplot(1, length(trialToPlot), j)
        errorbar_pn_hao(dataToplot);  % Plot the group average with error bars
        
        % Use a title that shows both trial name and block field
        title(sprintf('%s - %s', trialNames(trialToPlot(j)), blockField))
        ylim([0 1])
        ylabel('P(port entry)')
        xlim([0 size(dataToplot,2)])
        xticks([0 floor(size(dataToplot,2)/2) size(dataToplot,2)])
        xticklabels([win_time/-1 0 win_time])
        xlabel('Time from CS onset (s)')
    end
end







%% plot avg/group

color = {[0.8500, 0.3250, 0.0980];[0, 0.4470, 0.7410];[0.9290 0.6940 0.1250]};
% colors: 
% - color{1} = orange / shock
% - color{2} = blue / sucrose
% - color{3} = yellow / neutral

avg_fig = figure('Name','Average Port Entry Probability','NumberTitle','off');


% SUCROSE 

% get data
ty = 1;
binned_data = pe(ty).mxt.rew_ms;
avg_data = mean(binned_data);

% plot line
plot(avg_data)

% calculate error
err = std(binned_data)/sqrt(size(binned_data, 1));
neg_err = err * -1;
err_x = 1:size(avg_data,2);

% plot error
hold on
patch([err_x flip(err_x)], [avg_data-err flip(avg_data+err)], color{2}, 'FaceAlpha',0.25, 'EdgeColor','none')
hold off

% y-axis
ylabel('Average P(port entry)')
ylim([0 1])
yticks([0 0.5 1])

% x-axis
xlabel('Time from CS onset (s)')
xlim([0 size(avg_data,2)])
xticks([0 floor(size(avg_data,2)/2) size(avg_data,2)])
xticklabels([win_time/-1 0 win_time])

% legend
legend(trialNames(1), 'Location','northwest')


% NEUTRAL 
if isequal(trialToPlot,[1 4]) %rew/neut
    ty = 2;
elseif isequal(trialToPlot,[1 2 4]) %disc
    ty = 3;
else
    return;
end

% get data
binned_data = pe(ty).mxt.rew_ms;
avg_data = mean(binned_data);

% plot line
hold on
plot(avg_data, 'color', color{3})
hold off

% calculate error
err = std(binned_data)/sqrt(size(binned_data, 1));
neg_err = err * -1;
err_x = 1:size(avg_data,2);

% plot error
hold on
patch([err_x flip(err_x)], [avg_data-err flip(avg_data+err)], color{3}, 'FaceAlpha',0.25, 'EdgeColor','none')
hold off

% legend
legend(trialNames(1), '', trialNames(4), 'Location','northwest')


% SHOCK
if isequal(trialToPlot,[1 2 4])

    % get data
    ty = 2; % 
    binned_data = pe(ty).mxt.rew_ms;
    avg_data = mean(binned_data);
    
    % plot line
    hold on
    plot(avg_data, 'color', color{1})
    hold off
    
    % calculate error
    err = std(binned_data)/sqrt(size(binned_data, 1));
    neg_err = err * -1;
    err_x = 1:size(avg_data,2);
    
    % plot error
    hold on
    patch([err_x flip(err_x)], [avg_data-err flip(avg_data+err)], color{1}, 'FaceAlpha',0.25, 'EdgeColor','none')
    hold off
    
    % legend
    legend(trialNames(1), '', trialNames(4), '', trialNames(2), 'Location','northwest')
end

%% plot individual port latency

figure('Name','Port Latency (s)','NumberTitle','off')

for j = 1: length(trialToPlot) % looping for each CS
    ty = find(trialToPlot(j)==trialtypes);
    toPlot = pe(ty).latency.rew.txt;
    counter = 1;
    maxTrials = size(pe(1).latency.rew.txt, 1);

    for i = 1:mice % plotting latencies for 1 of the CS for all mice
        plotcoloum = ceil(mice/3);
        subplot(plotcoloum, ceil(mice/plotcoloum), counter)
        plot(toPlot(:, i))
        hold on
        title(AnimalIDcell(i,1))
        ylim([0 15])
        ylabel('Port Latency (s)')
        xlim([0 maxTrials])
        xticks([0 maxTrials/4 maxTrials/2 maxTrials*(3/4) maxTrials])
        xticklabels([0 maxTrials/4 maxTrials/2 maxTrials*(3/4) maxTrials])
        xlabel('Trial Number')
        legend(nm,'Location','northeast')
        counter = counter +1;
    end
end

%% plot avg/group port latency - run this then next section

% adding a mean matrix to the struct for each CS (grp mean per trial)
for i = 1: length(trialToPlot) % looping for each CS
    ty = find(trialToPlot(i)==trialtypes);
    rowMeans = mean(pe(ty).latency.rew.txt, 2);
    pe(ty).latency.rew.grp = rowMeans;
end

%% ^^ cont.

figure('Name','Average Port Latency (s)','NumberTitle','off')

for j = 1: length(trialToPlot) % looping for each CS
    ty = find(trialToPlot(j)==trialtypes);
    toPlot = pe(ty).latency.rew.grp;
    
    plot(toPlot)
    hold on
    ylabel('Average Port Latency (s)')
    ylim([0 15])
    xlabel('Trial Number')
    xlim([0 maxTrials])
    xticks([0 maxTrials/4 maxTrials/2 maxTrials*(3/4) maxTrials])
    xticklabels([0 maxTrials/4 maxTrials/2 maxTrials*(3/4) maxTrials])
    legend(nm, 'Location','northeast')
    hold on
end

%% plot avg port latency for each mouse 
% x-axis needs to be changed to 1-8 and make bar graph or dot graph

figure('Name','Average Port Latency (s) for each mouse','NumberTitle','off')


for j = 1: length(trialToPlot) % looping for each CS
    ty = find(trialToPlot(j)==trialtypes);
    toPlot = pe(ty).latency.rew.m;
    
    x = 1:mice;
    y = toPlot(1,:);
    scatter(x, y, 100, 'filled')

    hold on
    ylabel('Average Port Latency (s)')
    ylim([0 15])
    xlabel('Mouse ID')
    set(gca, 'XTick', 1:1:8)
    xlim([0 mice+1])
    xticklabels(AnimalIDcell(:,1))
    legend(nm, 'Location','northeast')
    hold on
end

%%  plot individual port latency - but w/ trials in sequential order

% create cell array with trial type and port latency for each mouse

% initialize 1x8 cell array
portLat = cell(1,8); 

% define first column values for each cell (trial type list)
firstColumnValues = trialtype{1,6};

% set total number of trials in the session
totalTrials = length(trialtype{1,6});

% fill each cell with a 120x2 matrix, setting the first column to the
% values above
for i = 1:mice
    portLat{i} = zeros(totalTrials,2);
    portLat{i}(:,1) = firstColumnValues;
end

% figure out second column values
for i = 1:mice
    sucCount = 1;
    neutCount = 1;
    shockCount = 1;
    for row = 1:totalTrials
        if portLat{i}(row,1) == 1
            portLat{i}(row,2) = pe(1).latency.rew.txt(sucCount,i);
            sucCount = sucCount + 1;
        elseif portLat{i}(row,1) == 2 & isequal(trialToPlot,[1 2 4])
            portLat{i}(row,2) = pe(2).latency.rew.txt(shockCount,i);
            shockCount = shockCount + 1;
        elseif portLat{i}(row,1) == 4 & isequal(trialToPlot,[1 2 4])
            portLat{i}(row,2) = pe(3).latency.rew.txt(neutCount,i);
            neutCount = neutCount + 1;
        elseif portLat{i}(row,1) == 4 & isequal(trialToPlot,[1 4])
            portLat{i}(row,2) = pe(2).latency.rew.txt(neutCount,i);
            neutCount = neutCount + 1;
        else
            error('Error: portLat could not be calculated properly.')
        end
    end
end

% check that this worked properly - on 2nd mouse as an example
if length(portLat{2}(:,2)) ~=  totalTrials
    error("Error: portLat check failed.")
end
%% 


% now onto the plotting

figure('Name','Port Latency (s)','NumberTitle','off')

xplot = 1:120;
counter = 1;
for i = 1:mice % looping for each mouse (for all CS) to plot line
    plotcoloum = ceil(mice/3);
    subplot(plotcoloum, ceil(mice/plotcoloum), counter)
    toPlot = portLat{i}(:,2);
    plot(toPlot, 'k')
    hold on
    title(AnimalIDcell(i,1))
    ylim([0 15])
    ylabel('Port Latency (s)')
    xlim([0 totalTrials])
    xticks([0 totalTrials/4 totalTrials/2 totalTrials*(3/4) totalTrials])
    xticklabels([0 totalTrials/4 totalTrials/2 totalTrials*(3/4) totalTrials])
    xlabel('Trial Number')
    counter = counter +1;
    hold on
    
    for row = 1:totalTrials
        if portLat{i}(row,1) == 1
            plotColor = color{2}; % blue dots for sucrose
        elseif portLat{i}(row,1) == 4
            plotColor = color{3}; % yellow dots for neutral
        elseif portLat{i}(row,1) == 2
            plotColor = color{1}; % red dots for shock
        end
    plot(xplot(row), portLat{i}(row,2), 'o', 'MarkerFaceColor', plotColor, 'MarkerEdgeColor', plotColor, 'MarkerSize', 3)
    end
end

% create liine objects for the legend
hSuc = line(nan, nan, 'Marker', 'o', 'Color', color{2}, 'MarkerFaceColor', color{2}, 'LineStyle', 'none', 'MarkerSize', 4.5); % Blue dot
hNeut = line(nan, nan, 'Marker', 'o', 'Color', color{3}, 'MarkerFaceColor', color{3}, 'LineStyle', 'none', 'MarkerSize', 4.5); % Yellow dot
hShock = line(nan, nan, 'Marker', 'o', 'Color', color{1}, 'MarkerFaceColor', color{1}, 'LineStyle', 'none', 'MarkerSize', 4.5); % Red dot


hold on
if isequal(trialToPlot,[1 4])
    legendHandle = legend([hSuc, hNeut], nm, 'Location', 'bestoutside');
else
    legendHandle = legend([hSuc, hShock, hNeut], nm, 'Location', 'bestoutside');
end
legendHandle.Position(1) = 0.6;
legendHandle.Position(2) = 0.28;

set(gcf, 'Position', [100, 200, 1800, 1100]);
%% 

%now try average



%% plot blocks (high rew vs high avers) for each CS type

color = {[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290 0.6940 0.1250]; [0.5, 0.8, 1]; [1.0, 0.7, 0.4]; [1.0, 0.85, 0.5]};
% colors: 
% - color{1} = blue / sucrose (or suc blks 1 + 2)
% - color{2} = orange / shock (or shk blks 1 + 2)
% - color{3} = yellow / neutral (or neut blks 1 + 2)
% - color{4} = light blue (suc blks 3 + 4)
% - color{5} = light orange (shk blks 3 + 4)
% - color{6} = light yellow (neut blks 3 + 4)

fig = figure;
nSubplots = 3; % bc 3 tones in disc

for j = 1:nSubplots % 1 suc, 2 shk, 3 neut
    subplot(1, 3, j);
    hold on;
    dataR = blkPE(j).blkR;
    dataA = blkPE(j).blkA;
    plot(dataR, 'LineStyle', '-', 'Color', color{j}); % high reward blocks avg
    plot(dataA, 'LineStyle', '--', 'Color', color{j}); % high aversive blocks avg
    title(trialNames(trialToPlot(j)));
    legend({'High Reward', 'High Averisve'}, 'Location', 'best')

    % x-axis
    xlim([0 size(dataA, 2)])
    xticks([0 floor(size(dataA,2)/2) size(dataA,2)])
    xticklabels([win_time/-1 0 win_time])
    xlabel('Time from CS onset (s)')

    % y-axis
    ylim([0 1])
    ylabel('P(port entry)')

    hold off;
end

%% plot blocks (all 4) for each CS type

fig2 = figure;
nSubplots = 3; % bc 3 tones in disc

for j = 1:nSubplots % 1 suc, 2 shk, 3 neut
    subplot(1, 3, j);
    hold on;
    data1 = blkPE(j).blk1;
    data2 = blkPE(j).blk2;
    data3 = blkPE(j).blk3;
    data4 = blkPE(j).blk4;
    plot(data1, 'LineStyle', '-', 'Color', color{j}); % blk 1 (high rew)
    plot(data2, 'LineStyle', '--', 'Color', color{j}); % blk 2 (high avers)
    plot(data3, 'LineStyle', '-', 'Color', color{j+3}); % blk 3 (high rew)
    plot(data4, 'LineStyle', '--', 'Color', color{j+3}); % blk 4 (high avers)
    title(trialNames(trialToPlot(j)));
    legend({'Block 1 (High Reward)', 'Block 2 (High Averisve)', 'Block 3 (High Reward)', 'Block 4 (High Aversive)'}, 'Location', 'best')

    % x-axis
    xlim([0 size(dataA, 2)])
    xticks([0 floor(size(dataA,2)/2) size(dataA,2)])
    xticklabels([win_time/-1 0 win_time])
    xlabel('Time from CS onset (s)')

    % y-axis
    ylim([0 1])
    ylabel('P(port entry)')

    hold off;
end

%% plot learned blocks (3+4 only) for each CS type

fig3 = figure;
nSubplots = 3; % bc 3 tones in disc

for j = 1:nSubplots % 1 suc, 2 shk, 3 neut
    subplot(1, 3, j);
    hold on;
    data3 = blkPE(j).blk3;
    data4 = blkPE(j).blk4;
    plot(data3, 'LineStyle', '-', 'Color', color{j}); % blk 3 (high rew)
    plot(data4, 'LineStyle', '--', 'Color', color{j}); % blk 4 (high avers)
    title(trialNames(trialToPlot(j)));
    legend({'Block 3 (High Reward)', 'Block 4 (High Aversive)'}, 'Location', 'best')

    % x-axis
    xlim([0 size(dataA, 2)])
    xticks([0 floor(size(dataA,2)/2) size(dataA,2)])
    xticklabels([win_time/-1 0 win_time])
    xlabel('Time from CS onset (s)')

    % y-axis
    ylim([0 1])
    ylabel('P(port entry)')

    hold off;
end

%% checking for sex diff 
% specifically for pavlovpilot2 > day12 > c1 = M, c2 = F

% create pe mouse x time matrix
for i = 1:length(trialToPlot)
    pe(i).mxtime = struct();
    data = pe(i).rawdata.rew; % get full data for one trial type
    data = data(:, 1501:3001, :); % take only second half (aka the 15s after tone begins)
    result = mean(data, 2); % get mean along 2nd dimension (time) -> making it 35x1x8
    result = squeeze(result); % making it 35x8
    pe(i).mxtime.rew = result;
end

% making table
for i = 1:length(trialToPlot)
    pe(i).mxtime.blks = zeros(4,mice);
end
    
% SUC blk 1
for j = 1:mice
    data = mean(pe(1).mxtime.rew(1:13, j));
    pe(1).mxtime.blks(1,j) = data;
end
% SUC blk 2
for j = 1:mice
    data = mean(pe(1).mxtime.rew(14:18, j));
    pe(1).mxtime.blks(2,j) = data;
end
% SUC blk 3
for j = 1:mice
    data = mean(pe(1).mxtime.rew(19:30, j));
    pe(1).mxtime.blks(3,j) = data;
end
% SUC blk 4
for j = 1:mice
    data = mean(pe(1).mxtime.rew(31:35, j));
    pe(1).mxtime.blks(4,j) = data;
end

% SHK blk 1
for j = 1:mice
    data = mean(pe(2).mxtime.rew(1:5, j));
    pe(2).mxtime.blks(1,j) = data;
end
% SHK blk 2
for j = 1:mice
    data = mean(pe(2).mxtime.rew(6:18, j));
    pe(2).mxtime.blks(2,j) = data;
end
% SHK blk 3
for j = 1:mice
    data = mean(pe(2).mxtime.rew(19:23, j));
    pe(2).mxtime.blks(3,j) = data;
end
% SHK blk 4
for j = 1:mice
    data = mean(pe(2).mxtime.rew(24:35, j));
    pe(2).mxtime.blks(4,j) = data;
end

% NEUT blk 1
for j = 1:mice
    data = mean(pe(3).mxtime.rew(1:5, j));
    pe(3).mxtime.blks(1,j) = data;
end
% NEUT blk 2
for j = 1:mice
    data = mean(pe(3).mxtime.rew(6:10, j));
    pe(3).mxtime.blks(2,j) = data;
end
% NEUT blk 3
for j = 1:mice
    data = mean(pe(3).mxtime.rew(11:15, j));
    pe(3).mxtime.blks(3,j) = data;
end
% NEUT blk 4
for j = 1:mice
    data = mean(pe(3).mxtime.rew(16:20, j));
    pe(3).mxtime.blks(4,j) = data;
end

% row 5 = high rew avg, row 6 = high avers avg
for i = 1:length(trialToPlot)
    for j = 1:mice
        avg1 = pe(i).mxtime.blks(1,j);
        avg3 = pe(i).mxtime.blks(3,j);
        pe(i).mxtime.blks(5,j) = (avg1 + avg3) / 2;
        avg2 = pe(i).mxtime.blks(2,j);
        avg4 = pe(i).mxtime.blks(4,j);
        pe(i).mxtime.blks(6,j) = (avg2 + avg4) / 2;
    end
end

%% plotting ^^
% for indiv mice

fig4 = figure;

for j = 1:nSubplots % 1 suc, 2 shk, 3 neut
    subplot(1, 3, j);
    hold on;
    x = 1:4;
    toPlot = pe(j).mxtime.blks(5,1:4); % males, high rew = blue, filled
    scatter(x, toPlot, 'b', 'filled'); 
    toPlot = pe(j).mxtime.blks(5,5:8); % females, high rew = pink, filled
    scatter(x, toPlot, 'm', 'filled');
    toPlot = pe(j).mxtime.blks(6,1:4); % males, high avers = blue, unfilled
    scatter(x, toPlot, 'b');
    toPlot = pe(j).mxtime.blks(6,5:8); % females, high avers = pink, unfilled
    scatter(x, toPlot, 'm');

    title(trialNames(trialToPlot(j)));
    legends = {'High Reward - Males', 'High Reward - Females', 'High Aversive - Males', 'High Aversive - Females'};
    legend(legends, 'Location', 'best')

    % x-axis
    xlim([1 4]);
    xticks(1:4);
    xlabel('Mouse Number');

    % y-axis
    ylim([0 1]);
    ylabel('P(port entry)');

    hold off;
end

sgtitle = 'Probabilty of Port Entry by State and Sex';

%% plotting ^^
% for avg of M vs F

fig5 = figure;

avgs = zeros(4, 3);
sems = zeros(4, 3);

% plotting scatter
for j = 1:length(trialToPlot) % 1 suc, 2 shk, 3 neut

    % male avg, high rew = blue, filled
    toPlot = pe(j).mxtime.blks(3,1:4); 
    avgs(1,j) = mean(toPlot);
    sems(1,j) = stderr(toPlot);
    scatter(j, mean(toPlot), 'b', 'filled'); 
    hold on;

    % female avg, high rew = pink, filled
    toPlot = pe(j).mxtime.blks(3,5:8); 
    avgs(2,j) = mean(toPlot);
    sems(2,j) = stderr(toPlot);
    scatter((j+0.2), mean(toPlot), 'm', 'filled');

    % male avg, high avers = blue, unfilled
    toPlot = pe(j).mxtime.blks(4,1:4); 
    avgs(3,j) = mean(toPlot);
    sems(3,j) = stderr(toPlot);
    scatter((j+0.4), mean(toPlot), 'b');

    % female avg, high avers = pink, unfilled
    toPlot = pe(j).mxtime.blks(4,5:8); 
    avgs(4,j) = mean(toPlot);
    sems(4,j) = stderr(toPlot);
    scatter((j+0.6), mean(toPlot), 'm');
end

% plotting error bars
errorbar(1:3, avgs(1,:), sems(1,:), 'LineStyle', 'none', 'Color', 'b', 'LineWidth', 1) % high rew - m
errorbar(1.2:3.2, avgs(2,:), sems(2,:), 'LineStyle', 'none', 'Color', 'm', 'LineWidth', 1) % high rew - f
errorbar(1.4:3.4, avgs(3,:), sems(3,:), 'LineStyle', 'none', 'Color', 'b') % high avers - m
errorbar(1.6:3.6, avgs(4,:), sems(4,:), 'LineStyle', 'none', 'Color', 'm') % high avers - f


title('Probabilty of Port Entry by State and Sex');
legend(legends, 'Location', 'best');

% x-axis
xlim([0.6 4]);
xticks(1.3:3.3);
xticklabels({'Sucrose CS', 'Shock CS', 'Neutral CS'})
xlabel('Trial Type')

% y-axis
ylim([0 1]);
ylabel('P(port entry)');


hold off;

%% plotting sex DIFFS

% create data array - row 7 = row 3 - row 4 (high rew - high avers (later blocks)) 
for j = 1:length(trialToPlot)
    for col = 1:mice
        pe(j).mxtime.blks(7,col) = pe(j).mxtime.blks(3,col) - pe(j).mxtime.blks(4,col);
    end
end

fig6 = figure;

avgs = zeros(2, 3);
sems = zeros(2, 3);

% plotting scatter
for j = 1:length(trialToPlot)

    % high rew - high avers : males (blue)
    toPlot = pe(j).mxtime.blks(7, 1:4);
    avgs(1,j) = mean(toPlot);
    sems(1,j) = stderr(toPlot);
    scatter(j, mean(toPlot), 'g', 'filled');
    hold on;

    % high rew - high avers : females (pink)
    toPlot = pe(j).mxtime.blks(7, 5:8);
    avgs(2,j) = mean(toPlot);
    sems(2,j) = stderr(toPlot);
    scatter((j+0.2), mean(toPlot), 'm', 'filled');
end

% plotting error bars 
errorbar(1:3, avgs(1,:), sems(1,:), 'LineStyle', 'none', 'Color', 'g', 'LineWidth', 1) % m
errorbar(1.2:3.2, avgs(2,:), sems(2,:), 'LineStyle', 'none', 'Color', 'm', 'LineWidth', 1) % f

title('Difference in P(PE) btw High Reward and High Aversive Blocks');
legends = {'Males', 'Females'};
legend(legends, 'Location', 'northeast');

% x-axis
xlim([0.6 3.6]);
xticks(1.1:3.1);
xticklabels({'Sucrose CS', 'Shock CS', 'Neutral CS'})
xlabel('Trial Type')

% y-axis
ylim([-0.4 0.4]);
yticks(-0.4:0.1:0.4)
ylabel('P(PE) in High Reward - P(PE) in High Aversive');

hold off;