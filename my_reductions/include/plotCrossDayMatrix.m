function plotCrossDayMatrix(matrix, dateList, titleStr, isPmatrix)
    % Plots cross-day decoding matrix or p-matrix as heatmap
    % matrix: nDays x nDays (accuracy or p-values)
    % dateList: cell array of dates (strings)
    % titleStr: string for title
    % isPmatrix: true/false → if true, plot p-values 0-1

    figure;
    if isPmatrix
        % Plot p-values 0-1
        imagesc(matrix, [0 1]);
        colormap(parula); % Good for p-values
        cbar = colorbar;
        cbar.Label.String = 'p-value';
    else
        % Plot accuracies 0-1
        imagesc(matrix, [0 1]);
        colormap(jet);
        cbar = colorbar;
        cbar.Label.String = 'Accuracy';
    end

    axis square;

    nDays = length(dateList);
    xticks(1:nDays);
    yticks(1:nDays);
    xticklabels(dateList);
    yticklabels(dateList);
    xtickangle(45);

    title(titleStr, 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Test Day');
    ylabel('Train Day');

    % Optional: overlay text
    for i = 1:nDays
        for j = 1:nDays
            if ~isnan(matrix(i,j))
                if isPmatrix
                    txt = sprintf('p=%.3f', matrix(i,j));
                    % Optionally highlight p<0.05
                    if matrix(i,j) < 0.05
                        text(j, i, txt, 'HorizontalAlignment','center','Color','r','FontSize',10,'FontWeight','bold');
                    else
                        text(j, i, txt, 'HorizontalAlignment','center','Color','w','FontSize',10);
                    end
                else
                    txt = sprintf('%.2f', matrix(i,j)*100);
                    text(j, i, txt, 'HorizontalAlignment','center','Color','w','FontSize',10);
                end
            end
        end
    end
end
