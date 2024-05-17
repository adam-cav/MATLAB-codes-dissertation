function [machNumbers, pulseSpeeds] = calcMachNumbers(posOnePeaksPos, posOneTime, posTwoPeaksPos, posTwoTime, distance)

signalCount = size(posOnePeaksPos, 1);

[machNumbers, pulseSpeeds] = deal(zeros(signalCount,1));

for signalIdx = 1:signalCount

    first = posOneTime(signalIdx,posOnePeaksPos(signalIdx,1));
    second = posTwoTime(signalIdx,posTwoPeaksPos(signalIdx,1));
    
    timeDifference = (second - first);
    pulseSpeeds(signalIdx) = distance / timeDifference;

    machNumbers(signalIdx) = pulseSpeeds(signalIdx) ./ 343;
    
end

% speed = distance / time