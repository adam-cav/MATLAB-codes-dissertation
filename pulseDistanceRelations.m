function [durationDistanceModel, pressureDistanceModel, modelDistanceVector, measuredDurations, measuredDistanceVector, modelPressures] = pulseDistanceRelations(data,time,distanceChange,peaksVals,peaksPos,troughsPos)

c0 = 343;
rho0 = 1.204;
gamma = 1.41;

distances = repmat(distanceChange, 1, ceil(size(peaksVals,2) / size(distanceChange,2)));
distances = distances(1:size(peaksVals,2)) / 1000;

distanceX = [0,distances];
distanceX = distanceX(1:size(peaksVals,2));

measuredDistanceVector = cumsum(distanceX);

modelDistanceVector = linspace(0,40,4000);

[durationDistanceModel, pressureDistanceModel] = deal(zeros(size(peaksVals,1),numel(modelDistanceVector)));
[measuredDurations, modelPressures] = deal(zeros(size(peaksVals,1),size(peaksVals,2)));

for signalIdx = 1:size(troughsPos,1)

    [p0, modelPressures(signalIdx,1)] = deal(max(data(signalIdx,:)));
    t0 = time(signalIdx,troughsPos(signalIdx,1)) - time(signalIdx,peaksPos(signalIdx,1));

    durationDistanceModel(signalIdx,:) = t0 .* sqrt(1+((gamma + 1)/2) .* ((p0 .* modelDistanceVector ./ (rho0 * c0^3 * t0))));
    pressureDistanceModel(signalIdx,:) = p0 ./ (sqrt(1+((gamma + 1)/2) .* ((p0 .* modelDistanceVector ./ (rho0 * c0^3 * t0)))));

    for peakIdx = 1:size(troughsPos,2)
        measuredDurations(signalIdx,peakIdx) = time(signalIdx,troughsPos(signalIdx,peakIdx)) - time(signalIdx,peaksPos(signalIdx,peakIdx));
    end

    for idx = 2:size(troughsPos,2)
        modelPressures(signalIdx,idx) = p0 ./ (sqrt(1+((gamma + 1)/2) .* ((p0 .* measuredDistanceVector(idx) ./ (rho0 * c0^3 * t0)))));
    end
end