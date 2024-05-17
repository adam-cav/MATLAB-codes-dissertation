function varargout = loadMultiAI2(commonFileName, filePath, transducerNumbers, signalLength, signalCount, time2pulse)

numTransducers = numel(transducerNumbers);
dataCell = cell(1, numTransducers);

for transIndex = 1:numTransducers
    % Initialise matrices for current transducer
    pressureValueArray = zeros(signalCount, signalLength + 1);
    timeValueArray = zeros(signalCount, signalLength + 1);
    
    transducer = transducerNumbers(transIndex);

    for signalIndex = 1:signalCount
        
        startIdx = 0;
        finishIdx = 0;
        
        % Generate file path    
        fileName = sprintf('%s%d',commonFileName,signalIndex);
        fileName = strrep(fileName, ' ', '');
        file = fullfile(filePath, fileName);

        % Load data
        data = readtable(file, 'VariableNamingRule', 'preserve');
        dataName = sprintf("AI %d/AI %d (Pa)", transducer, transducer);
        
        % Read pressure and time signals
        pressureVar = data.(dataName);
        timeVar = data.("Time (s)");

        % Trim signals
        riseSignal = diff(pressureVar);
        [~, maxElementPos] = max(riseSignal);
        startIdx = maxElementPos - time2pulse;
        finishIdx = startIdx + signalLength;
        
        pressureVar = pressureVar(startIdx:finishIdx);
        timeVar = timeVar(startIdx:finishIdx);

        % Store corresponding time array

        % Normalise time [0 -> N]
        % minTime = min(timeVar);
        % normalisedTimeVar = normalisedTimeVar - minTime;


        % Store data
        pressureValueArray(signalIndex, :) = pressureVar;
        timeValueArray(signalIndex, :) = timeVar;
        
    end
    
    % Store matrices for current transducer in cell 
    dataCell{transIndex} = struct('pressure', pressureValueArray, 'time', timeValueArray);
end

varargout = cell(1, (numel(dataCell) * 2)+1);

for i = 1:numel(dataCell)
    current = dataCell{i};

    pressureArray = current.pressure;
    timeArray = current.time;

    varargout{(i-1)*2 + 1} = pressureArray;
    varargout{(i-1)*2 + 2} = timeArray;
end
varargout{end} = 0:1/200000:(signalLength/200000);

% dataCell = splitCell(dataCell);