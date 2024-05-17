function [pressureDistanceModel, modelDistanceVector, measuredDistanceVector] = pulsePressureDistance(distanceChange,timeVector,peaksVals,peaksPos,troughsPos)

c0 = 343;
rho0 = 1.18;
gamma = 1.41;

distances = repmat(distanceChange, 1, ceil(size(peaksVals,2) / size(distanceChange,2)));
distances = distances(1:size(peaksVals,2)) / 1000;

distanceX = [0,distances];
distanceX = distanceX(1:size(peaksVals,2));

measuredDistanceVector = cumsum(distanceX);

modelDistanceVector = linspace(0,40,4000);

pressureDistanceModel = zeros(size(peaksVals,1),numel(modelDistanceVector));

for signalIdx = 1:size(peaksVals,1)
    p0 = peaksVals(signalIdx,1);
    t0 = timeVector(signalIdx,troughsPos(signalIdx,1)) - timeVector(signalIdx,peaksPos(signalIdx,1));

    pressureDistanceModel(signalIdx,:) = p0./(sqrt(1+((gamma + 1)/2) .* ((p0 .* modelDistanceVector ./ (rho0 * c0^3 * t0)))));
end