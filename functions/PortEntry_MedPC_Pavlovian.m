
clear all
close all

%% get files and port entries
[MedPCfile, folderDir] = GetFilesFromFolder(0,'.txt');
% MedPCfile = MedPCfile(contains(MedPCfile,'!'));
% folderDir = folderDir(contains(MedPCfile,'!'));

[allport, ~ , AnimalIDcell,trialtype,trialTS] = PortEntry_MedPC2mat(MedPCfile, folderDir);

%% get psth
win_time = 15;
win_half = win_time*100; % 10s before and after cue onset
binsize = 0.5; % 0.5s per bin
trialtypes = unique(trialtype{1}); % get trial types
for i = 1:length(trialtypes)
    ty = trialtypes(i);
    cue_onset = [trialtype{i}, trialTS{i}];
    pe(i) = processPE(cue_onset, win_half, binsize, allport, ty);
end
%save('port entry.mat','pe','trialNames', 'AnimalIDcell','allport');
%% plot individual 
session = 1;
mice = 3;
plotcoloum = ceil(mice/3);
trialNames = ["Sucrose CS"; "Shock CS"; "Neutral CS"];
trialToPlot = [1 3];
nm = trialNames(trialToPlot);
figure
for j = 1: length(trialToPlot)
    ty = find(trialToPlot(j)==trialtypes);
    dataToplot = pe(ty).mxt.rew_ms;
    counter = 1;
    for i = ((session-1)*mice)+1:(session)*mice
        subplot(plotcoloum, mice/plotcoloum, counter)
        plot(dataToplot(i,:,:))
        hold on
        title(AnimalIDcell(i,1))
        ylim([0 1])
        ylabel('P(port entry)')
        xlim([0 size(dataToplot,2)])
        xticks([0 floor(size(dataToplot,2)/2) size(dataToplot,2)])
        xticklabels([win_time/-1 0 win_time])
        xlabel('Time from CS onset')
        legend(nm,'Location','northwest')
        counter = counter +1;
    end
end

%% plot group
% figure
% 
% temp = ["Sucrose CS"; "Neutral CS"; "Shock CS"];
% label = ["No Stim"; "CS Stim"];
% color = {[0.8500, 0.3250, 0.0980];[0, 0.4470, 0.7410];[0.4940 0.1840 0.5560]};
% 
% for ty = 1:length(temp)
%     trialtypeTOplot = contains(trialNames,temp(ty));
%     trialName = trialNames(trialtypeTOplot);
%     trialtype = trialtypes(trialtypeTOplot);
%     subplot(1, 3, ty)
%     
%         dataToplot = pe(ty).mxt.rew_ms;
%         b(1) = errorbar_pn(1:size(dataToplot,2), mean(dataToplot), std(dataToplot)/sqrt(size(dataToplot, 1)), color{1});
%         hold on
%         dataToplot = pe(ty+3).mxt.rew_ms;
%         b(2) = errorbar_pn(1:size(dataToplot,2), mean(dataToplot), std(dataToplot)/sqrt(size(dataToplot, 1)), color{2});
%         ylim([0 1])
%         ylabel('P(port entry)')
%         xlim([0 size(dataToplot,2)])
%         xticks([0 floor(size(dataToplot,2)/2) size(dataToplot,2)])
%         xticklabels([-30 0 30])
%         xlabel('Time from CS onset')
%     
%     title(temp(ty))
%     legend(b,label)
% end











