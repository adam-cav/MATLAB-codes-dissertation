function [reflectionCoefficient, absorptionCoefficient] = energyReflAbsCoeff(data,peaksPos,troughsPos)
% Compute the complex frequency-dependent reflection coefficient of two pulses.
%
% reflectionCoefficient = energyReflCoeff(incidentPulse,reflectedPulse,lengthCorrection) 
% Loads the provided incident and reflected pulse, pads both arrays to desired 
% lengthCorrection, and returns the energy reflection coefficient.
%
% Inputs:
%  - incidentPulsePath: path/to/incident/pulse.mat
%  - reflectedPulsePath: path/to/reflected/pulse.mat
%  - lengthCorrection: desired length of pulse matrices
%
% Output:
%  - reflectionCoefficient: returns energy reflection coefficient.

reflectionCoefficient = zeros(size(data,1),size(peaksPos,2)-1);


for signalIdx = 1:size(data,1)
    for pulseIdx = 1:size(peaksPos,2)-1
        pulseOne = data(signalIdx,peaksPos(signalIdx,pulseIdx):troughsPos(signalIdx,pulseIdx));
        pulseTwo = data(signalIdx,peaksPos(signalIdx,pulseIdx+1):troughsPos(signalIdx,pulseIdx+1));

        reflectionCoefficient(signalIdx,pulseIdx) = sum(abs(pulseOne)).^2 / sum(abs(pulseTwo)).^2;
    end
end
absorptionCoefficient = 1-reflectionCoefficient;
