function pulseEnergyMatrix = pulseEnergies(data,peaksPos,troughsPos)
% Returns matrix of calculated pulse energies
signalCount = size(data,1);
numPeaks = size(peaksPos,2);

pulseEnergyMatrix = zeros(signalCount,numPeaks);

for signalIdx = 1:signalCount
    for pulseIdx = 1:numPeaks
        
        pulseEnergyMatrix(signalIdx,pulseIdx) = sum(abs(data(signalIdx,peaksPos(signalIdx,pulseIdx)-200:troughsPos(signalIdx,pulseIdx))).^2);

    end
end