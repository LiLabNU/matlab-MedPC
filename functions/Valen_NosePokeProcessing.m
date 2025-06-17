function [data, trialTS, Totals, AnimalIDcell, trial_average] = Valen_NosePokeProcessing(MedPCfile, folderDir, events, plot)

[data, trialTS, AnimalIDcell] = FRnosepoke_MedPC2mat(MedPCfile, folderDir);

%% raster plot for individual mice
bin = 60; %60s per line

if plot == "yes"
    imageRowCol = [ceil(size(MedPCfile,1)/8) 8]; %  8 rows for the subplot
    Plot_individualRaster_NP(data, trialTS, AnimalIDcell, imageRowCol, bin, events);
end
%% Compute # of events
%each row from Hao's outputs corresponds to one animal, in the list under
%'MedPCfile' to the right

PortEntry = data.PortEntry;
Sucrose = data.Sucrose;
ActiveNP = data.ActiveNP;
InactiveNP = data.InactiveNP;
for i=1:length(MedPCfile) % for i equal to 1 to the # of medPC files (i.e., animals)
    Totals.Sucrose(i,1)=sum(Sucrose(i,:));%set sucrose_n, in the i-th row and first column, to the sum of events in the i-th row of Sucrose
    Totals.ActiveNP(i,1)=sum(ActiveNP(i,:));
    Totals.InactiveNP(i,1)=sum(InactiveNP(i,:));
    Totals.PortEntry(i,1)=sum(PortEntry(i,:));
end


%
% if imageCol~=1 && imageRow~=1
%     Sucrose_n = reshape(Sucrose_n,[imageCol,imageRow]);
%     ActiveNP_n = reshape(ActiveNP_n,[imageCol,imageRow]);
%     InactiveNP_n = reshape(InactiveNP_n,[imageCol,imageRow]);
%     PortEntry_n = reshape(PortEntry_n,[imageCol,imageRow]);
% end
%% calculate and plot trial average
win_time = 30;
win_half = win_time*100; % 30s before and after cue onset
binsize = 0.5; % 0.5s per bin
CueTS = trialTS.CueTS;

for j = 1:length(events)
    event = cell2mat(events(j));
    dataToCalculate = data.(event);
    if CueTS{1} ~= 0
        trial_average.(event) = processPE(CueTS, win_half, binsize, dataToCalculate);
    else
        trial_average = [];
    end

    % plot individual trial average
    if plot == "plot"
        session = 1;
        mice = size(MedPCfile,1);
        plotcoloum = ceil(mice/4);

        dataToplot = trial_average.(event).mxt.rew_ms;
        counter = 1;
        for i = ((session-1)*mice)+1:(session)*mice
            subplot(plotcoloum, ceil(mice/plotcoloum), counter)
            plot(dataToplot(i,:,:))
            hold on
            title(AnimalIDcell(i,1))
            %ylim([-5 20])
            ylabel('Probability')
            xlim([0 size(dataToplot,2)])
            xticks([0 floor(size(dataToplot,2)/2) size(dataToplot,2)])
            xticklabels([win_time/-1 0 win_time])
            xlabel('Time from CS onset')
            if counter == 1
                legend(event,'Location','northwest')
            end
            counter = counter +1;
        end
    end
end
end
