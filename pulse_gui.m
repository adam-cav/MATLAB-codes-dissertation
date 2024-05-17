% function pulse_gui(data)
%     % Example matrix of signals
%     signals = data;  % 5 signals, 100 samples each
% 
%     % Create the main figure
%     fig = uifigure('Name', 'Signal and Pulse Selector', 'Position', [100, 100, 800, 600]);
% 
%     % Axes for plotting signals - Initialise before setting callbacks
%     ax = uiaxes(fig, 'Position', [180, 300, 600, 280]);
% 
%     % Dropdown for selecting signals - Initialise dropdown after axes
%     dd = uidropdown(fig, 'Position', [20, 550, 150, 30], 'Items', string(1:size(signals, 1)), ...
%                     'ValueChangedFcn', @(dd, event) plotSignal(dd, ax, signals, indices));
% 
%     % Table for start indices
%     startTable = uitable(fig, 'Position', [20, 20, 350, 250], 'ColumnName', {'Start Indices'});
% 
%     % Table for end indices
%     endTable = uitable(fig, 'Position', [430, 20, 350, 250], 'ColumnName', {'End Indices'});
% 
%     % Initialise indices structure properly to avoid scalar structure assignment error
%     indices = struct('starts', {cell(1, size(signals, 1))}, 'ends', {cell(1, size(signals, 1))});
%     for i = 1:size(signals, 1)
%         indices.starts{i} = [];
%         indices.ends{i} = [];
%     end
% 
%     % Plot the first signal by default
%     plotSignal(dd, ax, signals, indices);
% 
%     % Set up the click callback
%     ax.ButtonDownFcn = @(src, event) pickPoints(src, event, indices, startTable, endTable, signals, str2double(dd.Value));
% end
% 
% function plotSignal(dd, ax, signals, indices)
%     % Get the selected signal
%     idx = str2double(dd.Value);
% 
%     % Clear previous plots
%     cla(ax);
% 
%     % Plot the selected signal
%     plot(ax, signals(idx, :));
%     hold(ax, 'on');
%     if ~isempty(indices.starts{idx})
%         plot(ax, indices.starts{idx}, signals(idx, indices.starts{idx}), 'ro', 'MarkerFaceColor', 'r');
%     end
%     if ~isempty(indices.ends{idx})
%         plot(ax, indices.ends{idx}, signals(idx, indices.ends{idx}), 'bo', 'MarkerFaceColor', 'b');
%     end
%     hold(ax, 'off');
%     title(ax, ['Signal ' num2str(idx)]);
% end
% 
% function pickPoints(ax, event, indices, startTable, endTable, signals, idx)
%     % Get the point clicked
%     point = round(event.IntersectionPoint(1));
% 
%     % Determine the current state of the point
%     current_starts = indices.starts{idx};
%     current_ends = indices.ends{idx};
% 
%     % Check if the point is a start or end or needs to be removed
%     if ismember(point, current_starts)
%         current_starts(current_starts == point) = [];
%     elseif ismember(point, current_ends)
%         current_ends(current_ends == point) = [];
%     else
%         if isempty(current_starts) || (any(current_ends) && point < min(current_ends))
%             current_starts(end + 1) = point;  % Add as a new start
%         else
%             current_ends(end + 1) = point;  % Add as a new end
%         end
%     end
% 
%     % Sort indices after adding
%     current_starts = sort(current_starts);
%     current_ends = sort(current_ends);
% 
%     % Update indices
%     indices.starts{idx} = current_starts;
%     indices.ends{idx} = current_ends;
% 
%     % Update tables
%     startTable.Data = indices.starts{idx};
%     endTable.Data = indices.ends{idx};
% 
%     % Redraw the plot
%     plotSignal(dd, ax, signals, indices);
% end
% function pulse_gui(data)
%     % Example matrix of signals
%     signals = data;  % data contains signals, each row is a signal
% 
%     % Create the main figure
%     fig = uifigure('Name', 'Signal and Pulse Selector', 'Position', [100, 100, 800, 600]);
% 
%     % Axes for plotting signals - Initialise before setting callbacks
%     ax = uiaxes(fig, 'Position', [180, 300, 600, 280]);
% 
%     % Dropdown for selecting signals
%     dd = uidropdown(fig, 'Position', [20, 550, 150, 30], 'Items', string(1:size(signals, 1)), ...
%                     'ValueChangedFcn', @(dd, event) plotSignal(dd, ax, signals, indices));
% 
%     % Initialise mode buttons
%     btnStart = uibutton(fig, 'push', 'Text', 'Select Start Points', ...
%                         'Position', [20, 500, 150, 30], 'ButtonPushedFcn', @(btn,event) setMode('start'));
%     btnEnd = uibutton(fig, 'push', 'Text', 'Select End Points', ...
%                       'Position', [20, 460, 150, 30], 'ButtonPushedFcn', @(btn,event) setMode('end'));
% 
%     % Initialise indices structure properly
%     indices = struct('starts', {cell(1, size(signals, 1))}, 'ends', {cell(1, size(signals, 1))});
%     for i = 1:size(signals, 1)
%         indices.starts{i} = [];
%         indices.ends{i} = [];
%     end
% 
%     % Set initial mode
%     mode = 'start';  % default mode
% 
%     % Define mode-setting function
%     function setMode(newMode)
%         mode = newMode;
%     end
% 
%     % Table for start indices
%     startTable = uitable(fig, 'Position', [20, 20, 350, 250], 'ColumnName', {'Start Indices'});
% 
%     % Table for end indices
%     endTable = uitable(fig, 'Position', [430, 20, 350, 250], 'ColumnName', {'End Indices'});
% 
%     % Plot the first signal by default
%     plotSignal(dd, ax, signals, indices);
% 
%     % Set up the click callback
%     ax.ButtonDownFcn = @(src, event) pickPoints(src, event, indices, startTable, endTable, signals, str2double(dd.Value), mode);
% end
% 
% function plotSignal(dd, ax, signals, indices)
%     % Get the selected signal
%     idx = str2double(dd.Value);
% 
%     % Clear previous plots
%     cla(ax);
% 
%     % Plot the selected signal
%     plot(ax, signals(idx, :));
%     hold(ax, 'on');
%     if ~isempty(indices.starts{idx})
%         plot(ax, indices.starts{idx}, signals(idx, indices.starts{idx}), 'ro', 'MarkerFaceColor', 'r');
%     end
%     if ~isempty(indices.ends{idx})
%         plot(ax, indices.ends{idx}, signals(idx, indices.ends{idx}), 'bo', 'MarkerFaceColor', 'b');
%     end
%     hold(ax, 'off');
%     title(ax, ['Signal ' num2str(idx)]);
% end
% 
% function pickPoints(ax, event, indices, startTable, endTable, signals, idx, mode)
%     % Get the point clicked
%     point = round(event.IntersectionPoint(1));
% 
%     % Determine the current state of the point and mode
%     if strcmp(mode, 'start')
%         modifyIndices(indices.starts{idx}, point);
%     elseif strcmp(mode, 'end')
%         modifyIndices(indices.ends{idx}, point);
%     end
% 
%     function modifyIndices(indexList, point)
%         if ismember(point, indexList)
%             indexList(indexList == point) = [];  % Remove the point
%         else
%             indexList(end + 1) = point;  % Add the point
%         end
%         indexList = sort(indexList);  % Sort after modifying
%     end
% 
%     % Update indices and tables
%     indices.starts{idx} = sort(indices.starts{idx});
%     indices.ends{idx} = sort(indices.ends{idx});
%     startTable.Data = indices.starts{idx};
%     endTable.Data = indices.ends{idx};
% 
%     % Redraw the plot
%     plotSignal(dd, ax, signals, indices);
% end
