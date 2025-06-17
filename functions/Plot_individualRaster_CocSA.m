function Plot_individualRaster_CocSA(data, AnimalIDcell, imageRowCol, bin, legendLabels)

A = data.PumpOn;
B = data.NoCocActiveNP;
C = data.InactiveNP;
D = data.CocTrigNP;
%E = data.PortEntry;
totalmice = size(A,1) ;


figure('Units', 'normalized', 'Position', [0,0,1,1]); % Create a full-screen figure

for i = 1:totalmice
    ax = subplot(imageRowCol(1), imageRowCol(2), i); % Get the handle of the subplot


    a = reshape(A(i,:), [size(A(i,:),2)/bin, bin])';
    b = reshape(B(i,:), [size(B(i,:),2)/bin, bin])';
    c = reshape(C(i,:), [size(C(i,:),2)/bin, bin])';
    d = reshape(D(i,:), [size(D(i,:),2)/bin, bin])';

    % Collect handles from plotRaster_hao
    h1 = plotRaster_hao(a, "dot", 'k.', 10, ax);
    h2 = plotRaster_hao(b, "dot", 'g.', 10, ax);
    h3 = plotRaster_hao(c, "dot", 'c.', 10, ax);
    h4 = plotRaster_hao(d, "dot", 'b.', 10, ax);

    hold on;
    if i == 1
        % Dummy plots for legend
        hLegend1 = plot(nan, nan, 'k.', 'MarkerSize', 10); % Black dot for 'Reward'
        hLegend2 = plot(nan, nan, 'g.', 'MarkerSize', 10); % Green dot for 'ActiveNP'
        hLegend3 = plot(nan, nan, 'c.', 'MarkerSize', 10); % Green dot for 'InactiveNP'
        hLegend4 = plot(nan, nan, 'b.', 'MarkerSize', 10);  % Red line for 'RewardTrigNP'

        legend([hLegend1, hLegend2, hLegend3, hLegend4], legendLabels, 'Location', 'northeast');

        % Hide the dummy plots
        set([hLegend1, hLegend2, hLegend3, hLegend4], 'Visible', 'off');

        xlabel('Time (ms)');
        ylabel('Time (min)');
    end
    hold off;
    xlim([0 size(a,2)])
    ylim([0 size(a,1)+1])
    title(AnimalIDcell(i,1) + ' ' + AnimalIDcell(i,2) + ' ' + AnimalIDcell(i,3));
end

end
