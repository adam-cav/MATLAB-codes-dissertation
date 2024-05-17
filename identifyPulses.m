function [pks, pklocs, trs, trlocs] = identifyPulses(data, minPeakDistance, checkSize)
%IDENTIFYPULSES Detect peaks and troughs in signals.
% Inputs:
%   data - Matrix where each row represents a signal.
%   minPeakDistance - Minimum distance between detected peaks.
%   checkSize - Half-width of the window for local minimum detection.

% Get the number of signals
sigCount = size(data,1);

% Determine peak properties from the first signal to establish array sizes
temp = findpeaks(data(1,:),'MinPeakProminence',0.12*max(data(1,:)), 'MinPeakDistance',minPeakDistance);
peakCount = numel(temp);

% Initialise arrays for peaks and troughs
[pks, pklocs, tempPks] = deal(zeros(sigCount, peakCount)); % [values, locations]
[trs, trlocs] = deal(zeros(sigCount, peakCount-1)); % [values, locations]

for i = 1:sigCount
    % Detect peaks with specified prominence and distance
    % [pks(i,:), pklocs(i,:)] = findpeaks(data(i,:),'MinPeakProminence',0.12*max(data(i,:)), 'MinPeakDistance',minPeakDistance, 'NPeaks', peakCount);
    [pksTemp, pklocsTemp] = findpeaks(data(i,:),'MinPeakProminence',0.12*max(data(i,:)), 'MinPeakDistance',minPeakDistance, 'NPeaks', peakCount);
    if numel(pksTemp) ~= numel(size(pks,2))
        count = size(pksTemp);
        pks(i,1:count(end)) = pksTemp;
        pklocs(i,1:count(end)) = pklocsTemp;
    else
        [pks(i,:), pklocs(i,:)] = findpeaks(data(i,:),'MinPeakProminence',0.12*max(data(i,:)), 'MinPeakDistance',minPeakDistance, 'NPeaks', peakCount);
    end

    % Detect pulse start points based on max rise in signal
    for j = 1:max(count)
        checkWin = data(i,pklocs(i,j)-checkSize:pklocs(i,j)+checkSize);
        [~, maxRiseIdx] = max(diff(checkWin));
        trueMaxRiseIdx = maxRiseIdx + (pklocs(i,j) - checkSize);
        tempPks(i,j) = pklocs(i,j);
        pklocs(i,j) = trueMaxRiseIdx;
    end

    for j = 1:peakCount - 1
        % Define search start immediately after the current peak
        start = tempPks(i,j) + 1;
        % Determine the search end point
        if j < peakCount
            finish = tempPks(i,j+1);  
        else
            finish = numel(data(i,:));  
        end

        foundTrough = false;
        for ij = start:min(finish, numel(data(i,:)) - checkSize)
            if (ij + checkSize) > numel(data(i,:)) || (ij - checkSize) < 1
                continue; % Skip if the window exceeds data bounds
            end
            localWindow = data(i, ij - checkSize:ij + checkSize);
            [localMin, localIdx] = min(localWindow);

            if localIdx == checkSize + 1
                trlocs(i,j) = ij;
                trs(i,j) = localMin;
                foundTrough = true;
                break; % Exit after finding the first local minimum
            end
        end

        if ~foundTrough && (finish > start)
            % Fallback to find minimum in the entire range if no local trough is detected
            [trs(i,j), idx] = min(data(i, start:finish));
            trlocs(i,j) = start + idx - 1; % Adjust index due to 1-based indexing
        end
    end
end
end



% function [pks, pklocs, trs, trlocs] = identifyPulses(data,minPeakDistance,checkSize)
% 
% sigCount = size(data,1);
% 
% temp = findpeaks(data(1,:),'MinPeakProminence',0.12*max(data(1,:)), 'MinPeakDistance',minPeakDistance);
% 
% [pks, pklocs, trs, trlocs] = deal(zeros(sigCount,size(temp,2)));
% 
% for i = 1:sigCount
% 
%     [pks(i,:), pklocs(i,:)] = findpeaks(data(i,:),'MinPeakProminence',0.12*max(data(i,:)), 'MinPeakDistance',minPeakDistance);
% 
%     for j = 1:numel(pks(i,:)) - 1
% 
%         start = pklocs(i,j) + 1;
% 
%         if j < numel(pklocs)
%             finish = pklocs(j+1);  
%         else
%             finish = numel(data(i,:));  
%         end
% 
%         foundTrough = false;
%         for ij = start:min(finish, numel(data(i,:)) - checkSize)
%             localWindow = data(i, ij - checkSize:ij + checkSize);
%             [localMin, localIdx] = min(localWindow);
% 
%             if localIdx == checkSize + 1
% 
%                 trlocs(i,j) = ij;
%                 trs(i,j) = localMin;
%                 foundTrough = true;
%                 break
% 
%             end
%         end
% 
%         if ~foundTrough && (finish > start)
% 
%             [trs(i,j),index] = min(data(i,start:finish));
%             trlocs(i,j) = start + index - 1;
% 
%         end
%     end
% 
% end