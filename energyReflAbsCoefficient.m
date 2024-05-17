function [energyReflCoeff, energyAbsCoeff] = energyReflAbsCoefficient(data,pulseEnergyMatrix,peaksPos)
% Returns the energy reflection and absorption coefficients

signalCount = size(data,1);
numPeaks = size(peaksPos,2);

[energyReflCoeff, energyAbsCoeff] = deal(zeros(signalCount,numPeaks-1));

for i = 1:signalCount
    for j = 1:numPeaks-1
        energyReflCoeff(i,j) = pulseEnergyMatrix(i,j+1) / pulseEnergyMatrix(i,j);
        energyAbsCoeff(i,j) = abs(1 - energyReflCoeff(i,j));
    end
end