function [reflectionCoefficient, frequencyVector, incidentFFT, reflectionCoefficientMagnitude, reflectionCoefficientPhase, timeVector] = freqReflCoeff(incidentPulsePath, reflectedPulsePath, lengthCorrection, Fs, numSignals)
% Compute the complex frequency-dependent reflection coefficient of two pulses.
%
% [reflCoeff,~,~,~,~] = freqReflCoeff(incidentPulse,reflectedPulse,lengthCorrection,Fs) 
% Loads the provided incident and reflected pulse, pads both arrays to the desired 
% lengthCorrection, and returns the complex frequency-dependent reflection coefficient.
%
% Inputs:
%  - incidentPulsePath: path/to/incident/pulse.mat
%  - reflectedPulsePath: path/to/reflected/pulse.mat
%  - lengthCorrection: desired length of pulse matrices
%  - Fs: sampling frequency for frequency and time vectors 
%
% Output:
%  - reflectionCoefficient: returns complex reflection coefficient for each frequency in frequencyVector
%
% Optional outputs:
%  - frequencyVector = returns the frequency vector
%  - timeVector = returns the time vector 
%  - reflectionCoefficientMagnitude: returns abs(reflectionCoefficient)  
%  - reflectionCoefficientPhase: returns angle(reflectionCoefficient)  

[reflectionCoefficient,reflectionCoefficientMagnitude, reflectionCoefficientPhase] = deal(zeros(numSignals,lengthCorrection/2));
incidentArray = zeros(1,lengthCorrection);
reflectedArray = zeros(1,lengthCorrection);

basePath = regexp(incidentPulsePath,'\d','once');
str1 = incidentPulsePath(1:basePath-1);
str2 = incidentPulsePath(basePath+1:end);
str3 = reflectedPulsePath(basePath+1:end);

for i = 1:numSignals
    incidentFileName = sprintf('%s%d%s',str1,i,str2);
    reflectedFileName = sprintf('%s%d%s',str1,i,str3);

    T = 1/Fs;
    
    incidentFile = load(incidentFileName);
    incidentName = fieldnames(incidentFile);
    incidentArray = incidentFile.(incidentName{1});
    incidentArray = [incidentArray, zeros(1,lengthCorrection-numel(incidentArray))];
    
    reflectedFile = load(reflectedFileName);
    reflectedName = fieldnames(reflectedFile);
    reflectedArray = reflectedFile.(reflectedName{1});
    reflectedArray = [reflectedArray, zeros(1,lengthCorrection-numel(reflectedArray))];
    
    L = numel(incidentArray);
    
    incidentFFT = fft(incidentArray);
    inc2 = abs(incidentFFT/L);
    inc1 = inc2(1:L/2+1);
    inc1(2:end-1) = 2*inc1(2:end-1);
    incidentFFT = inc1(1:end-1);

    reflectedFFT = fft(reflectedArray);
    refl2 = abs(reflectedFFT/L);
    refl1 = refl2(1:L/2+1);
    refl1(2:end-1) = 2*refl1(2:end-1);
    reflectedFFT = refl1(1:end-1);

    timeVector = (0:L-1)*T; 
    frequencyVector = Fs/L*(0:L-1); frequencyVector = frequencyVector(1:L/2);
    
    reflectionCoefficient(i,:) = reflectedFFT.^2 ./ incidentFFT.^2;
    reflectionCoefficientMagnitude(i,:) = abs(reflectionCoefficient(i,:));
    reflectionCoefficientPhase(i,:) = angle(reflectionCoefficient(i,:));

end