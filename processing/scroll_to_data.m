function scroll_to_data(x, col7, col8)
    % Create a figure with a plot
    fig = figure;
    ax = axes(fig);
    plot(ax, x, '-o'); % '-o' makes the points clickable
    title(ax, 'Click a point to see its data');

    % Enable data cursor mode
    dcm = datacursormode(fig);
    set(dcm, 'UpdateFcn', @(~, event_obj) displayData(event_obj));
    set(dcm, 'Enable', 'on', 'SnapToDataVertex', 'on'); % Require click

    function txt = displayData(event_obj)
        % Get the index of the clicked point
        idx = event_obj.DataIndex;

        % Determine the range of data to display around the clicked point
        range = max(1, idx):min(idx+30, length(x));

        % Display the data in a new figure
        dataFig = findobj('Type', 'figure', 'Name', 'Data Display');
        if isempty(dataFig)
            dataFig = figure('Name', 'Data Display', 'NumberTitle', 'off');
        else
            figure(dataFig);
        end
        clf(dataFig);
        figPos = get(dataFig, 'Position');
        dataTbl = uitable(dataFig, 'Data', [(range)', x(range), col7(range), col8(range)], ...
                          'ColumnName', {'Index', 'X', 'Col7', 'Col8'}, ...
                          'ColumnEditable', [false, true, true, true], ...
                          'Position', [20, 20, figPos(3)-40, figPos(4)-40], ...
                          'CellEditCallback', @(src, event) updateData(src, event, range));

        % Return text to display in the data cursor (optional)
        txt = {['Index: ', num2str(idx)], ...
               ['X: ', num2str(x(idx))]};
    end

    function updateData(src, event, range)
        % Update the original data array with the edited value
        rowIndex = range(event.Indices(1));
        if event.Indices(2) == 2 % X column
            x(rowIndex) = event.NewData;
        elseif event.Indices(2) == 3 % Col7 column
            col7(rowIndex) = event.NewData;
        elseif event.Indices(2) == 4 % Col8 column
            col8(rowIndex) = event.NewData;
        end
        assignin('base', 'updated_col_1', x);
        assignin('base', 'updated_col_2', col7);
        assignin('base', 'updated_col_3', col8);
    end
end
