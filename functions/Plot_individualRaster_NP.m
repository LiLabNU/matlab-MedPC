function Plot_individualRaster_NP(data, trialTS, AnimalIDcell, imageRowCol, bin, legendLabels)

PortEntry = data.PortEntry;
Sucrose = data.Sucrose;
ActiveNP = data.ActiveNP;
InactiveNP = data.InactiveNP;
totalmice = size(PortEntry,1) ;

trialTS = trialTS.Cue;


figure('Units', 'normalized', 'Position', [0,0,1,1]); % Create a full-screen figure

for i = 1:totalmice
    ax = subplot(imageRowCol(1), imageRowCol(2), i); % Get the handle of the subplot
    if isempty(trialTS{i})

        pe = reshape(PortEntry(i,:), [size(PortEntry(i,:),2)/bin, bin])';
        np = reshape(ActiveNP(i,:), [size(ActiveNP(i,:),2)/bin, bin])';
        inactnp = reshape(InactiveNP(i,:), [size(InactiveNP(i,:),2)/bin, bin])';
        suc = reshape(Sucrose(i,:), [size(Sucrose(i,:),2)/bin, bin])';

    else
        win_time = 30;
        win_half = win_time*100; % 30s before and after cue onset
        for j = 1: size(trialTS{i},1)-1
            currentTrial = trialTS{i}(j);
            temp = PortEntry(i,:);
            pe(j,:) = temp(currentTrial-win_half:currentTrial+win_half);
            temp = ActiveNP(i,:);
            np(j,:) = temp(currentTrial-win_half:currentTrial+win_half);
            temp = InactiveNP(i,:);
            inactnp(j,:) = temp(currentTrial-win_half:currentTrial+win_half);
            temp = Sucrose(i,:);
            suc(j,:) = temp(currentTrial-win_half:currentTrial+win_half);
        end
    end
    % Collect handles from plotRaster_hao
    h1 = plotRaster_hao(pe, "dot", 'k.', 8, ax);
    h2 = plotRaster_hao(np, "dot", 'g.', 8, ax);
    h3 = plotRaster_hao(inactnp, "dot", 'c.', 8, ax);
    h4 = plotRaster_hao(suc, "line", 'r', 4, ax);

    hold on;
    if i == 1
        % Dummy plots for legend
        hLegend1 = plot(nan, nan, 'k.', 'MarkerSize', 10); % Black dot for 'Port Entry'
        hLegend2 = plot(nan, nan, 'g.', 'MarkerSize', 10); % Green dot for 'Active Nose Poke'
        hLegend3 = plot(nan, nan, 'c.', 'MarkerSize', 10); % Green dot for 'Inactive Nose Poke'
        hLegend4 = plot(nan, nan, 'r.', 'MarkerSize', 10);  % Red line for 'Sucrose'

        legend([hLegend1, hLegend2, hLegend3, hLegend4], legendLabels, 'Location', 'northeast');

        % Hide the dummy plots
        set([hLegend1, hLegend2, hLegend3, hLegend4], 'Visible', 'off');

        xlabel('Time (ms)');
        ylabel('Time (min)');
    end
    hold off;
    xlim([0 size(pe,2)])
    ylim([0 size(pe,1)+1])
    title(AnimalIDcell(i,1) + ' ' + AnimalIDcell(i,2) + ' ' + AnimalIDcell(i,3));
end

end
