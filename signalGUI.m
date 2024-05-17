function signalGUI(data)
    % Create figure and UI components
    fig = figure('Name', 'Signal Picker GUI', 'Position', [100, 100, 800, 600]);
    btn_select_signal = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Select Signal', 'Position', [20, 550, 150, 30], 'Callback', @selectSignal);
    ax = axes('Parent', fig, 'Position', [0.1, 0.1, 0.8, 0.7]);
    txt_selected_points = uicontrol(fig, 'Style', 'edit', 'Position', [20, 20, 300, 200], 'Max', 2, 'Min', 0, 'HorizontalAlignment', 'left');
    
    % Global variables
    selected_signal = [];
    clicked_points = [];
    hover_text = [];
    
    % Callback functions
    function selectSignal(~, ~)
        if ~isempty(data)
            % Implement signal selection logic
            % For simplicity, let's assume data is a matrix where each row represents a signal
            prompt = 'Enter the index of the signal to select:';
            dlgtitle = 'Select Signal';
            dims = [1 35];
            definput = {'1'};
            answer = inputdlg(prompt, dlgtitle, dims, definput);
            if ~isempty(answer)
                index = str2double(answer{1});
                if ~isnan(index) && index >= 1 && index <= size(data, 1)
                    selected_signal = data(index, :);
                    plotSignal();
                else
                    msgbox('Invalid signal index.', 'Error', 'error');
                end
            end
        else
            msgbox('No data available.', 'Error', 'error');
        end
    end
    
    function plotSignal
        plot(ax, selected_signal);
        xlabel('Time');
        ylabel('Amplitude');
        title('Selected Signal');
        % Set up callback for mouse clicks on the plot
        set(ax, 'ButtonDownFcn', @plotClickCallback);
    end

    function plotClickCallback(~, event)
        % Get coordinates of the clicked point
        x = event.IntersectionPoint(1);
        y = event.IntersectionPoint(2);
        % Store the clicked point
        clicked_points = [clicked_points; x, y];
        % Update text box with selected points
        updateSelectedPoints();
        % Plot clicked points
        if isempty(hover_text)
            hold(ax, 'on');
            hover_text = text(x, y, ['(', num2str(x, '%.2f'), ', ', num2str(y, '%.2f'), ')'], 'Color', 'red', 'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
            hold(ax, 'off');
        else
            set(hover_text, 'Position', [x, y]);
            set(hover_text, 'String', ['(', num2str(x, '%.2f'), ', ', num2str(y, '%.2f'), ')']);
        end
    end

    function updateSelectedPoints
        % Update text box with selected points
        if ~isempty(clicked_points)
            points_str = cellstr(num2str(clicked_points, '(%.2f, %.2f)\n'));
            set(txt_selected_points, 'String', points_str);
        else
            set(txt_selected_points, 'String', '');
        end
    end
end
