function varargout = splitCell(cellStructure)
    varargout = cell(1, numel(cellStructure) * 2);

    for i = 1:numel(cellStructure)
        currentStruct = cellStructure{i};

        pressureArray = currentStruct.pressure;
        timeArray = currentStruct.time;

        varargout{(i-1)*2 + 1} = pressureArray;
        varargout{(i-1)*2 + 2} = timeArray;
    end
end