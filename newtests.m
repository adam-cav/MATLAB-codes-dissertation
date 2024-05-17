clc
clear 
close all

%% Constants

gamma  = 1.401; % Specific heat ratio at 20C
P0 = 1.013*10^5; % atmospheric pressure at 20C [Pa]
p = linspace(1,80000,79999); % Pressure array [Pa]
gasConstR = 287.05; % Specific gas constant of air [J/kg K]
temperature = 18 + 273.15; % Air temperature [k]
rho0 = P0 / (gasConstR * temperature); % density of air
c0 = 343; % sound speed at 20C [m/s]
c = c0 + ((gamma - 1) / 2) * (c0 / P0) .* p; % pressure dependent sound speed [m/s]

%% loading

Fs = 200000;
path = '/Users/adamcavanagh/Desktop/Uni/Year 3/Project/Data/exports_090424/txt/';
transducers = 2:4;
path2 = './txt/';
transducers2 = 1;
transducers3 = [1,3];

signalLength = round(0.111 * Fs); % Signal length after loading
signalCount = 5; % Number of ruptures 
samples2firstPulse = 750; % number of samples till first pulse

name_e40 = 'e40_000'; % 40 um mylar
name_e23 = 'e23_000'; % 23 um mylar
name_ePa = 'ePa_000'; % single layer paper

name_emptyTube = 'emptyTube_000';
name_data = 'mylar_40_';

name_a6_40 = 'a6-40_000';
name_a6_23 = 'a6-23_000';
name_a6_Pa = 'a6-Pa_000';

name_a12_40 = 'a12-40_000';
name_a12_23 = 'a12-23_000';
name_a12_Pa = 'a12-Pa_000';

name_a18_40 = 'a18-40_000';
name_a18_23 = 'a18-23_000';
name_a18_Pa = 'a18-Pa_000';

name_a25_40 = 'a25-40_000';
name_a25_23 = 'a25-23_000';
name_a25_Pa = 'a25-Pa_000';

[e_40p1, e_40t1, e_40p2, e_40t2, e_40p3, e_40t3, normalTime] = loadMultiAI2(name_e40,path,transducers,signalLength,signalCount,samples2firstPulse);
[e_23p1, e_23t1, e_23p2, e_23t2, e_23p3, e_23t3] = loadMultiAI2(name_e23,path,transducers,signalLength,signalCount,samples2firstPulse);
[e_Pap1, e_Pat1, e_Pap2, e_Pat2, e_Pap3, e_Pat3] = loadMultiAI2(name_ePa,path,transducers,signalLength,signalCount,samples2firstPulse);

[emptyTubep, emptyTubet] = loadMultiAI2(name_emptyTube,path2,transducers2,signalLength,signalCount,samples2firstPulse);
[datatap1, datatat1, normTimeData] = loadMultiAI2(name_data,path2,1,18000,5,samples2firstPulse);

[a6_40p1, a6_40t1, a6_40p2, a6_40t2, a6_40p3, a6_40t3] = loadMultiAI2(name_a6_40,path,transducers,signalLength,signalCount,samples2firstPulse);
[a6_23p1, a6_23t1, a6_23p2, a6_23t2, a6_23p3, a6_23t3] = loadMultiAI2(name_a6_23,path,transducers,signalLength,signalCount,samples2firstPulse);
[a6_Pap1, a6_Pat1, a6_Pap2, a6_Pat2, a6_Pap3, a6_Pat3] = loadMultiAI2(name_a6_Pa,path,transducers,signalLength,signalCount,samples2firstPulse);

[a12_40p1, a12_40t1, a12_40p2, a12_40t2, a12_40p3, a12_40t3] = loadMultiAI2(name_a12_40,path,transducers,signalLength,signalCount,samples2firstPulse);
[a12_23p1, a12_23t1, a12_23p2, a12_23t2, a12_23p3, a12_23t3] = loadMultiAI2(name_a12_23,path,transducers,signalLength,signalCount,samples2firstPulse);
[a12_Pap1, a12_Pat1, a12_Pap2, a12_Pat2, a12_Pap3, a12_Pat3] = loadMultiAI2(name_a12_Pa,path,transducers,signalLength,signalCount,samples2firstPulse);

[a18_40p1, a18_40t1, a18_40p2, a18_40t2, a18_40p3, a18_40t3] = loadMultiAI2(name_a18_40,path,transducers,signalLength,signalCount,samples2firstPulse);
[a18_23p1, a18_23t1, a18_23p2, a18_23t2, a18_23p3, a18_23t3] = loadMultiAI2(name_a18_23,path,transducers,signalLength,signalCount,samples2firstPulse);
[a18_Pap1, a18_Pat1, a18_Pap2, a18_Pat2, a18_Pap3, a18_Pat3] = loadMultiAI2(name_a18_Pa,path,transducers,signalLength,signalCount,samples2firstPulse);

[a25_40p1, a25_40t1, a25_40p2, a25_40t2, a25_40p3, a25_40t3] = loadMultiAI2(name_a25_40,path,transducers,signalLength,signalCount,samples2firstPulse);
[a25_23p1, a25_23t1, a25_23p2, a25_23t2, a25_23p3, a25_23t3] = loadMultiAI2(name_a25_23,path,transducers,signalLength,signalCount,samples2firstPulse);
[a25_Pap1, a25_Pat1, a25_Pap2, a25_Pat2, a25_Pap3, a25_Pat3] = loadMultiAI2(name_a25_Pa,path,transducers,signalLength,signalCount,samples2firstPulse);

%% 
% instance = 'e'; % e(mpty), a6, a12, a18, a25
% amplitude = 'Pa'; % Pa(per), 23 (mylar), 40 (mylar)
% micPos = 3; % mic position
% sig = 2;
numPeaks = [];
% distanceChange = [6800, 3800]; % pos.1
% distanceChange = [6280, 4320]; % pos.2
distanceChange = [3550, 7050]; % pos.3

% data = eval(sprintf('%s_%sp%d',instance,amplitude,micPos));
% time = eval(sprintf('%s_%st%d',instance,amplitude,micPos));
% data = datatap1;
% time = datatat1;

% % Find minimum between peaks
% for i = 1:numel(peakLocs)-1
%     [troughVals(i), troughIndex] = min(data(sig,peakLocs(i):peakLocs(i+1)));
%     troughLocs(i) = troughIndex + peakLocs(i) - 1;
% end
% 
% searchWin = min(length(data(sig,:)),peakLocs(end) + minPeakDistance);
% [troughVals(end), troughIndex] = min(data(sig,peakLocs(end):searchWin));
% troughLocs(end) = troughIndex + peakLocs(end) - 1;

% % Find local minimum after peak
% [peakVals, peakLocs, peakw, peakp] = findpeaks(data(sig,:),'MinPeakProminence',0.12*max(data(sig,:)), 'MinPeakDistance', minPeakDistance);
% % [troughVals, troughLocs,troughw,troughp] = findpeaks(-data(sig,:),'MinPeakProminence',0.1*max(-data(sig,:)), 'MinPeakDistance', minPeakDistance);
% [troughVals, troughLocs] = deal(zeros(size(peakLocs)));
% 
% windowSize = 500;
% for j = 1:numel(peakLocs) - 1
% 
%     start = peakLocs(j) + 1;
% 
%     if j < numel(peakLocs)
%         finish = peakLocs(j+1);  % Correctly assign finish to the next peak's location
%     else
%         finish = numel(data(sig,:));  % If it's the last peak, finish at the end of the array
%     end
% 
%     foundTrough = false;
%     for ij = start:min(finish, numel(data(sig,:)) - windowSize)
% 
%         localWindow = data(sig, ij - windowSize:ij+windowSize);
%         [localMin, localIdx] = min(localWindow);
% 
%         if localIdx == windowSize + 1
%             troughLocs(j) = ij;
%             troughVals(j) = localMin;
%             foundTrough = true;
%             break
%         end
% 
%     end
% 
%     if ~foundTrough && (finish > start)
% 
%         [troughVals(j),idx] = min(data(sig,start:finish));
%         troughLocs(j) = start + idx - 1;
% 
%     end
% 
% end

minPeakDistance = 1800; % samples
checkSize = 400; % samples

% [peakVals, peakLocs, troughVals, troughLocs] = identifyPulses(data,minPeakDistance,checkSize);
% [peakVals, peakLocs, troughVals, troughLocs] = extractPeaks(data,time,7,minPeakDistance,500,distanceChange);

% energy = pulseEnergies(data,peakLocs(:,1:end-1),troughLocs);
% avgEnergy = mean(energy,1)
% boxData = [energy(:,1)', avgEnergy(1); energy(:,2)', avgEnergy(2)];
% boxData = [energy(:,1)'; energy(:,2)'];



%% 

crossSectionArea = pi * 0.0255^2;
energyTrapz = zeros(5,5);
for i = 1:5
    for j = 1:5
        energyTrapz(i,j) = trapz(1/Fs,rho0 * crossSectionArea*(data(i,emptyPeakLocs(i,j):emptyTroughLocs(i,j+1)).^2));
    end
end 
%% 
[~, ~, ~, durations] = pulseDistanceRelations(data,time,distanceChange,emptyPeakVals,emptyPeakLocs,emptyTroughLocs);
energy = pulseEnergies(data,emptyPeakLocs(:,1:end-1),emptyTroughLocs);
[reflCoeff, absCoeff] = energyReflAbsCoefficient(data,energy,emptyPeakLocs(:,1:end-1));


[stanDevA, meansA] = std(emptyPeakVals);
% [stanDevB, meansB] = std(durations);
[stanDevC, meansC] = std(reflCoeff);
[stanDevD, meansD] = std(absCoeff);
% [stanDevE, meansE] = std(energyTrapz);
[stanDevF, meansF] = std(energy);

clc

disp(['Average peaks: ', sprintf(' & %1.0f',round(meansA(1:5)))])
% disp(['Energies: ', num2str(round(meansE))])
% disp(['Integrated energies: ', sprintf(' & %0.2f',meansE(1:5))])
disp(['Energies: ', sprintf(' & %0.3e',meansF(1:5))])
disp(['STD: ', sprintf(' & %0.3e',stanDevF)])
disp(['Refl: ', sprintf(' & %0.3f',meansC(1:5))])
disp(['STD: ', sprintf(' & %0.3f',stanDevC)])
disp(['Abs: ', sprintf(' & %0.3f', meansD(1:5))])
disp(['STD: ', sprintf(' & %0.3f', stanDevD)])
% disp(['Standard deviation: ', num2str(stanDevA)])
% disp(['Incident duration: ', num2str(meansB(1)*1000)])
% disp(['Standard deviation: ', num2str(stanDevB(1)*1000)])

%% FFTs

adjustment = 0;
sig = 5;
instance = 'a6';
instance1 = 'e';
amplitude = '40';
micPos = 3;
emptyData = eval(sprintf('%s_%sp%d',instance1,amplitude,micPos));
data = eval(sprintf('%s_%sp%d',instance,amplitude,micPos));
% time = eval(sprintf('%s_%st%d',instance,amplitude,micPos));

[emptyPeakVals, emptyPeakLocs, emptyTroughVals, emptyTroughLocs] = identifyPulses(emptyData,minPeakDistance,checkSize);
[dataPeakVals, dataPeakLocs, dataTroughVals, dataTroughLocs] = identifyPulses(data,minPeakDistance,checkSize);
% [peakVals2, peakLocs2, troughVals2, troughLocs2] = identifyPulses(a18_Pap3,minPeakDistance,checkSize);
% [peakVals3, peakLocs3, troughVals3, troughLocs3] = identifyPulses(a12_Pap3,minPeakDistance,checkSize);
% [peakVals4, peakLocs4, troughVals4, troughLocs4] = identifyPulses(a6_Pap3,minPeakDistance,checkSize);

% [a6_40peaksVals, a6_40peaksPos, a6_40troughsVals, a6_40troughsPos] = identifyPulses(a6_40p3,minPeakDistance,checkSize);
% [a12_40peaksVals, a12_40peaksPos, a12_40troughsVals, a12_40troughsPos] = identifyPulses(a12_40p3,minPeakDistance,checkSize);
% [a18_40peaksVals, a18_40peaksPos, a18_40troughsVals, a18_40troughsPos] = identifyPulses(a18_40p3,minPeakDistance,checkSize);
% [a25_40peaksVals, a25_40peaksPos, a25_40troughsVals, a25_40troughsPos] = identifyPulses(a25_40p3,minPeakDistance,checkSize);
% 
% mylar40peaks = [a6_40peaksVals(:,1:7); a12_40peaksVals(:,1:7); a18_40peaksVals(:,1:7); a25_40peaksVals];
% [mylar40peakStd, mylar40MeanPeaks] = std(mylar40peaks); 

energyData = pulseEnergies(data,emptyPeakLocs(:,1:end-1),emptyTroughLocs);
[reflCoeffData, absCoeffData] = energyReflAbsCoefficient(data,energyData,emptyPeakLocs(:,1:end-1));

energyEmpty = pulseEnergies(emptyData,emptyPeakLocs(:,1:end-1),emptyTroughLocs);
[reflCoeffEmpty, absCoeffEmpty] = energyReflAbsCoefficient(emptyData,energyEmpty,dataPeakLocs(:,1:end-1));

emptyEnergyChange = energyEmpty(sig,1) - energyEmpty(sig,2);
dataEnergyChange = energyData(sig,1) - energyData(sig,2);

meanEnergyEmpty = mean(energyEmpty,1);
meanEnergyData = mean(energyData,1);

avgEmptyEnergyChange = meanEnergyData(1) - meanEnergyData(2);
avgDataEnergyChange = meanEnergyEmpty(1) - meanEnergyEmpty(2);

[stanDevC, meansC] = std(reflCoeffData);
[stanDevD, meansD] = std(absCoeffData);
[stanDevF, meansF] = std(energyData);

emptyPeakValsdB = 20*log10(emptyPeakVals/2e-5);
dataPeakValsdB = 20*log10(dataPeakVals/2e-5);

dataMeans = mean(dataPeakVals,1);
dataMeansdB = 20*log10(mean(dataPeakVals,1)/2e-5);
emptyMeans = mean(emptyPeakVals,1);
emptyMeansdB = 20*log10(mean(emptyPeakVals,1)/2e-5);

emptyPressureRefl = emptyMeans(2) / emptyMeans(1);
emptyPressureAbs = abs(1 - emptyPressureRefl);

dataPressureRefl = dataMeans(2) / dataMeans(1);
dataPressureSampleRefl = dataPressureRefl / emptyPressureRefl;

reflCoeffEmpty = mean(reflCoeffEmpty,1);
absCoeffEmpty = mean(absCoeffEmpty,1);
reflCoeffData = mean(reflCoeffData,1);


energyRefl = reflCoeffData(1) / reflCoeffEmpty(1);

for i = 1:5
    emptyPressureRefl(i) = emptyMeans(i+1) / emptyMeans(i);
end

clc

% disp(['ABH peak: ', sprintf('%0.3f ',peakVals1(sig,1:2))])
% disp(['ABH peak dB: ', sprintf('%0.3f ',peakVals1dB(sig,1:2))])
% disp(['ABH mean peak: ', sprintf('%0.3f ', means(1:2))])
% disp(['In dB: ', sprintf('%0.3f ', meansdB(1:2))])
% 
% disp(['ABH p2p: ', num2str(abs((peakVals(sig,1)-peakVals(sig,2)) - (peakVals1(sig,1)-peakVals1(sig,2))))])
% disp(['ABH p2p dB: ', num2str(abs((peakValsdB(sig,1) - peakValsdB(sig,2)) - (peakVals1dB(sig,1) - peakVals1dB(sig,2))))])
% disp(['ABH p2p avg: ', num2str((means(1) - means(2)) - (means1(1) - means1(2)))])
% disp(['ABH p2p avg dB: ', num2str((meansdB(1) - meansdB(2)) - (means1dB(1) - means1dB(2)))])
% disp(['ABH refl: ', num2str(reflCoeffData(1))])
% disp(['ABH abs: ', num2str(absCoeffData(1))])
% disp(['Refl: ', sprintf('%0.3f', abs(energyRefl))])
% disp(['Abs: ', sprintf('%0.3f', abs(energyAbs))])
% disp(['Energy reduction: ', sprintf('%1.3e',abs(emptyEnergyChange-dataEnergyChange))])
% disp(['Avg energy reduction: ', sprintf('%1.3e',abs(avgEmptyEnergyChange-avgDataEnergyChange))])
% disp(['ABH pressure R: ', sprintf('%1.3f',dataPressureRefl)])
% disp(['ABH pressure a: ', sprintf('%1.3f',dataPressureAbs)])

% disp(['Empty peak: ', sprintf('%0.3f ',peakVals(sig,1:2))])
% disp(['Empty peak dB: ', sprintf('%0.3f ',peakValsdB(sig,1:2))])
% disp(['Empty energies: ' sprintf('& %1.3e ',meanEnergyEmpty(1:5))])
% disp(['ABH mean peak: ', sprintf('%1.0f ', dataMeans(1:2))])
% disp(['In dB: ', sprintf('%1.1f ', dataMeansdB(1:2))])
disp(['ABH p2p avg: ',  sprintf('%1.0f ',(dataMeans(1) - dataMeans(2)) - (emptyMeans(1) - emptyMeans(2)))])
disp(['ABH p2p avg dB: ',  sprintf('%1.3f ',(dataMeansdB(1) - dataMeansdB(2)) - (emptyMeansdB(1) - emptyMeansdB(2)))])

sprintf('%s_%sp%d',instance1,amplitude,micPos)
sprintf('%s_%sp%d',instance,amplitude,micPos)
disp(['Empty energy R: ', sprintf('& %1.3f ',reflCoeffEmpty(1:5))])
% disp(['Empty energy R: ', sprintf('%1.3f ',absCoeffEmpty(1))])
disp(['Empty pressure R: ', sprintf('& %1.3f ',emptyPressureRefl)])
% disp(['Empty pressure a: ', sprintf('%1.3f ',emptyPressureAbs)])
disp(['ABH energy R: ', sprintf('%0.3f ', reflCoeffData(1))])
disp(['ABH true energy R: ', sprintf('%0.3f ', energyRefl)])
% disp(['ABH true energy a: ', sprintf('%0.3f ', energyAbs)])
disp(['ABH pressure R: ', sprintf('%1.3f ',dataPressureRefl)])
disp(['ABH true pressure R: ', sprintf('%1.3f ',dataPressureSampleRefl)])
% disp(['ABH true pressure a: ', sprintf('%1.3f ',dataPressureAbs)])
% disp(['Empty mean peak: ', sprintf('%1.0f ', emptyMeans(1:2))])
% disp(['ABH mean peak: ', sprintf('%1.0f ', dataMeans(1:2))])
% disp(['In dB: ', sprintf('%1.1f ', means1dB(1:2))])
% disp(['Empty abs: ', sprintf('%1.3f',absCoeffEmpty(1))])
% disp(['Empty pressure R: ', sprintf('%1.3f',emptyPressureRefl)])
% disp(['Empty pressure a: ', sprintf('%1.3f',emptyPressureAbs)])

return

% [~,NFFTs1,~,~,freq] = pulseFFT(a25_40p3,2,peakLocs1,troughLocs1,Fs);
% [~,NFFTs2] = pulseFFT(a18_40p3,2,peakLocs2,troughLocs2,Fs);
% [~,NFFTs3] = pulseFFT(a12_40p3,2,peakLocs3,troughLocs3,Fs);
% [~,NFFTs4] = pulseFFT(a6_40p3,2,peakLocs4,troughLocs4,Fs);

pulse = 1;
offset = 3; % for pulse plots, in samples 
length = 1500;

% [value1, index1] = max(wdenoise(data(sig,peakLocs1(sig,pulse)-offset:end)));
% [value2, index2] = max(wdenoise(data(sig,peakLocs1(sig,pulse+1)-offset:end)));
% [value3, index3] = max(wdenoise(emptyData(sig,peakLocs(sig,pulse)-offset:end)));
% [value4, index4] = max(wdenoise(emptyData(sig,peakLocs(sig,pulse+1)-offset:end)));

close all
% figure('Position',[2, 500, 800, 650])
figure
tiledlayout(1,1,"Padding",'tight')
nexttile 
hold on 
% h1 = plot(wdenoise(data(sig,peakLocs1(sig,pulse)-offset:end)),'LineWidth',0.5,'Color',[162 20 47]/255);
h2 = plot(wdenoise(data(sig,dataPeakLocs(sig,pulse+1)-offset:end)),'LineWidth',0.5,'Color','r');
% scatter(index1,value1,40,'sb','LineWidth',1)
% scatter(index2,value2,40,'sb','LineWidth',1)
h3 = plot(wdenoise(emptyData(sig,emptyPeakLocs(sig,pulse)-offset:end)),'LineWidth',0.5,'Color','b');
h4 = plot(wdenoise(emptyData(sig,emptyPeakLocs(sig,pulse+1)-offset:end)),'LineWidth',0.5,'Color',[0.4660 0.6740 0.1880]);
% scatter(index3,value3,40,'sk','LineWidth',1)
% scatter(index4,value4,40,'sk','LineWidth',1)
xlabel('Time, ms')
ylabel('Pressure, Pa')
xticks(0:100:1500)
xticklabels((0:100:1500) * 1000 / Fs)
xlim([0 1200])
box on 
grid on
legend([h3 h4 h2],'Empty incident','Empty reflected','Sample one reflected')
% text(index1, value1, sprintf(' %1.0f Pa',value1),"Position",[index1+7, value1],'Color','b','FontSize',11)
% text(index2, value2, sprintf(' %1.0f Pa',value2),"Position",[index2+7, value2+500],'Color','b','FontSize',11)
% text(index2, value2, [sprintf('%1.0f Pa ',value2) '\rightarrow' ],"Position",[index2-10, value2],'Color','b','HorizontalAlignment','right','FontSize',10)
% text(index3, value3, sprintf(' %1.0f Pa',value3),"Position",[index3+7, value3],'FontSize',11)
% text(index4, value4, sprintf(' %1.0f Pa',value4),"Position",[index4+7, value4+180],'FontSize',11)
set(gca,'FontSize',17)


return
figure
TL = tiledlayout(2,1,"Padding",'tight');
nexttile
hold on
plot(normalTime,wdenoise(data(sig,:),14) - adjustment*mean(data(sig,:)),'-b')
xlim([0 0.1])
scatter(normalTime(dataPeakLocs(sig,1:end-2)),dataPeakVals(sig,1:end-2) - adjustment*mean(data(sig,:)),75,'_r','LineWidth',0.75)
scatter(normalTime(dataTroughLocs(sig,1:end-1)),dataTroughVals(sig,1:end-1) - adjustment*mean(data(sig,:)),75,'_r','LineWidth',0.75)
xlabel('Time, s')
ylabel('Pressure, Pa')
sigAx = gca;
sigAx.FontSize = 15;
grid on
box on

nexttile
hold on
plot(normalTime,wdenoise(emptyData(sig,:),14) - adjustment*mean(emptyData(sig,:)),'-b')
xlim([0 0.1])
scatter(normalTime(emptyPeakLocs(sig,1:end-2)),emptyPeakVals(sig,1:end-2) - adjustment*mean(emptyData(sig,:)),75,'_r','LineWidth',0.75)
scatter(normalTime(emptyTroughLocs(sig,1:end-1)),emptyTroughVals(sig,1:end-1) - adjustment*mean(emptyData(sig,:)),75,'_r','LineWidth',0.75)
xlabel('Time, s')
ylabel('Pressure, Pa')
sigAx = gca;
sigAx.FontSize = 15;
grid on
box on

%% 
[~, dataPeakLocs, ~, ~] = identifyPulses(e_40p2,minPeakDistance,checkSize);
[~, peakLocs2, ~, ~] = identifyPulses(e_40p3,minPeakDistance,checkSize);

[mNums, pulseSpeeds] = calcMachNumbers(dataPeakLocs,e_40t2,peakLocs2,e_40t3,1.37)

%% 

paperRollOff = spectralRolloffPoint(paperFFT(1:2000)',freq(1:2000));
m23RollOff = spectralRolloffPoint(mylar23FFT(1:2000)',freq(1:2000));
m40RollOff = spectralRolloffPoint(mylar40FFT(1:2000)',freq(1:2000));

figure
tiledlayout(1,1,'Padding','tight')
nexttile
hold on
h1 = plot(freq,paperFFT,'-b',LineWidth=1);
xline(paperRollOff,'-r',[num2str(paperRollOff), ' Hz'],'LineWidth',0.5,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','right')
h2 = plot(freq,mylar23FFT,'--b',LineWidth=1.3);
xline(m23RollOff,'--r',[num2str(m23RollOff), ' Hz'],'LineWidth',0.5,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','right')
h3 = plot(freq,mylar40FFT,':b',LineWidth=2);
xline(m40RollOff,':r',[num2str(m40RollOff), ' Hz'],'LineWidth',0.5,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','left')
xlim([0 2000])
ylim([0 1.0025])
xlabel('Frequency, Hz')
ylabel('Amplitude   ^{|Y|}/_{max(|Y|)}')
legend([h1 h2 h3], 'Paper', '23 \mum','40 \mum')
set(gca,'FontSize',17)
box on

return
close all

% Plot signal with peaks and troughs
figure(Position=[400 400 660 700])
TL = tiledlayout(5,1,"Padding",'tight');
for i = 1:5
    nexttile
    hold on
    plot(normalTime,wdenoise(data(i,:),14),'-b')
    xlim([0 0.1])
    % ylim([0 15500])
    scatter(normalTime(emptyPeakLocs(i,1:end-2)),emptyPeakVals(i,1:end-2),75,'_r','LineWidth',0.75)
    scatter(normalTime(emptyTroughLocs(i,1:end-1)),emptyTroughVals(i,1:end-1),75,'_r','LineWidth',0.75)
    sigAx = gca;
    sigAx.FontSize = 15;
    grid on
    box on 
end
xlabel(TL,'Time, s','FontSize',15)
ylabel(TL,'Pressure, Pa','FontSize',15)

figure
TL = tiledlayout(1,1,"Padding",'tight');
nexttile
hold on
plot(normalTime,wdenoise(data(sig,:),14) - adjustment*mean(data(sig,:)),'-b')
xlim([0 0.1])
scatter(normalTime(emptyPeakLocs(sig,1:end-2)),emptyPeakVals(sig,1:end-2) - adjustment*mean(data(sig,:)),75,'_r','LineWidth',0.75)
scatter(normalTime(emptyTroughLocs(sig,1:end-1)),emptyTroughVals(sig,1:end-1) - adjustment*mean(data(sig,:)),75,'_r','LineWidth',0.75)
xlabel('Time, s')
ylabel('Pressure, Pa')
sigAx = gca;
sigAx.FontSize = 15;
grid on
box on

%% 
amplitudes = {'Pa';'23';'40'};
titles = {'100 GSM paper'; '23 \mum Mylar'; '40 \mum Mylar'};

% ylims = [-1000 20000; -5000 40000; -10000 60000];
% figure(Position=[500 50 1200 350])
% signalTL = tiledlayout(1,3,"TileSpacing",'tight','Padding','tight');

figure(Position=[500 50 1000 350])
energyTL = tiledlayout(1,3,"TileSpacing",'tight','Padding','tight');

for i = 1:3
    % if i == 1
    %     data = eval(sprintf('%s_%sp%d',instance, string(amplitudes(i)), micPos));
    %     for j = 1:5
    %         data(j,:) = data(j,:) - 0.5*mean(data(j,:));
    %     end
    %     time = eval(sprintf('%s_%st%d',instance, string(amplitudes(i)), micPos));
    % else
        data = eval(sprintf('%s_%sp%d',instance, string(amplitudes(i)), micPos));
        time = eval(sprintf('%s_%st%d',instance, string(amplitudes(i)), micPos));
    % end
    [emptyPeakVals, emptyPeakLocs, emptyTroughVals, emptyTroughLocs] = identifyPulses(data,minPeakDistance,checkSize);
    
    nexttile(signalTL)
    hold on
    plot(normalTime,wdenoise(data(sig,:)),'-b')
    xlim([0 0.1])
    ylim(ylims(i,:))
    scatter(normalTime(emptyPeakLocs(sig,1:end-1)),emptyPeakVals(sig,1:end-1),'xr','LineWidth',2)
    scatter(normalTime(emptyTroughLocs(sig,1:end-1)),emptyTroughVals(sig,1:end-1),'xr','LineWidth',2)
    title(titles(i),'FontWeight','normal','Interpreter','tex'); %,Interpreter='latex')
    sigAx = gca;
    sigAx.FontSize = 15;
    grid on
    hold off

    % energy = pulseEnergies(data,peakLocs(:,1:end-1),troughLocs);
    % boxData = [energy(:,1)'; energy(:,2)'];
    % nexttile(energyTL)
    % b1 = bar(boxData');
    % b1(1).FaceColor = 'k';
    % b1(2).FaceColor = [192, 192, 192] ./ 255;
    % title(titles(i),'FontWeight','demi','Interpreter','tex'); %,Interpreter='latex')
    % enrgAx = gca;
    % enrgAx.FontSize = 15;
end

% xlabel(signalTL,'Time, s','FontSize',15)
% ylabel(signalTL,'Pressure, Pa','FontSize',15)

xlabel(energyTL,'Rupture','FontSize',15)
ylabel(energyTL,'Pulse Energy, \Sigma |P_i|^2',Interpreter='tex',fontsize=15)

%% 
[dataPeakVals, dataPeakLocs, ~, dataTroughLocs] = identifyPulses(e_Pap3,minPeakDistance,checkSize);
[peakVals2, peakLocs2, ~, troughLocs2] = identifyPulses(e_23p3,minPeakDistance,checkSize);
[peakVals3, peakLocs3, ~, troughLocs3] = identifyPulses(e_40p3,minPeakDistance,checkSize);

[durationModel1, pressureModel1, modelDistance, measuredDurations1, distanceVector] = pulseDistanceRelations(e_Pap3,e_Pat3,distanceChange,dataPeakVals,dataPeakLocs,dataTroughLocs);
[durationModel2, pressureModel2, ~, measuredDurations2] = pulseDistanceRelations(e_23p3,e_23t3,distanceChange,peakVals2,peakLocs2,troughLocs2);
[durationModel3, pressureModel3, ~, measuredDurations3] = pulseDistanceRelations(e_40p3,e_40t3,distanceChange,peakVals3,peakLocs3,troughLocs3);


xtikMax = ceil(distanceVector(1:7)/10)*10;
majorTicks = 0:2.5:max(xtikMax);

close all
edistanceFig = figure;
% edistanceFig.Position = [902,3,1544,404];
eTL = tiledlayout(1,2,"TileSpacing",'tight','Padding','tight');
nexttile 
hold on 
h1 = plot(modelDistance,mean(durationModel1,1),'-r');
h2 = plot(modelDistance,mean(durationModel2,1),'--r');
h3 = plot(modelDistance,mean(durationModel3,1),':r');
h4 = scatter(distanceVector,mean(measuredDurations1,1),40,'rs','filled');
h5 = scatter(distanceVector,mean(measuredDurations2,1),40,'r^','filled');
h6 = scatter(distanceVector,mean(measuredDurations3,1),40,'rv','filled');
xticks(majorTicks)
xlim([0 35])
grid on 
box on
legend([h1 h2 h3], 'Paper', '23 \mum','40 \mum',Location='northwest')
ylabel('Pulse duration, s','FontSize',16)
eAx1 = gca;

nexttile
hold on
h1 = plot(modelDistance,mean(pressureModel1,1),'-b');
h2 = plot(modelDistance,mean(pressureModel2,1),'--b');
h3 = plot(modelDistance,mean(pressureModel3,1),':b');
h4 = scatter(distanceVector,mean(dataPeakVals,1),'bs','filled');
h5 = scatter(distanceVector,mean(peakVals2,1),'b^','filled');
h6 = scatter(distanceVector,mean(peakVals3,1),'bv','filled');
xticks(majorTicks)
xlim([0 35])
grid on 
box on
legend([h1 h2 h3], 'Paper', '23 \mum','40 \mum',Location='northeast')
xlabel(eTL,'Distance travelled, m','FontSize',16)
ylabel('Pressure, Pa','FontSize',16)
eAx2 = gca;
set([eAx1 eAx2], 'FontSize',16)
%% 

[dataPeakVals, ~, ~, ~] = identifyPulses(a25_40p3,minPeakDistance,checkSize);
[peakVals2, ~, ~, ~] = identifyPulses(a18_40p3,minPeakDistance,checkSize);
[peakVals3, ~, ~, ~] = identifyPulses(a12_40p3,minPeakDistance,checkSize);
[peakVals4, ~, ~, ~] = identifyPulses(a12_40p3,minPeakDistance,checkSize);
[peakVals5, ~, ~, ~] = identifyPulses(e_40p3,minPeakDistance,checkSize);

[firstSTD, firstMean] = std(dataPeakVals,[],1);
[secondSTD, secondMean] = std(peakVals2,[],1);
[thirdSTD, thirdMean] = std(peakVals3,[],1);
[fourSTD, fourMean] = std(peakVals4,[],1);
[fiveSTD, fiveMean] = std(peakVals5,[],1);

%%

close all 

sig = 3;


figure(position=[725 400 525 420])
tiledlayout(1,1,Padding='tight')
nexttile
hold on
h1 = errorbar(firstMean(1:5)/1000,firstSTD(1:5)/1000, 'Color','b','LineWidth',1);
h2 = errorbar(secondMean(1:5)/1000,secondSTD(1:5)/1000, 'Color','r','LineWidth',1);
% h3 = errorbar(thirdMean(1:5)/1000,thirdSTD(1:5)/1000, 'Color',[233 134 64] / 255);
% h4 = errorbar(fourMean(1:5)/1000,fourSTD(1:5)/1000, 'Color',[124 154 92] / 255);
h5 = errorbar(fiveMean(1:5)/1000,fiveSTD(1:5)/1000, 'Color',[165 0 165] / 255,'LineWidth',1);
% xticklabels({'i';'r_1';'r_2';'r_3';'r_4';'r_5';'r_6'})
xticks([1 2 3 4 5 6])
xticklabels({'i';'r_1';'r_2';'r_3';'r_4'})
xlim([0.95 5.05])
xlabel('Pulse Instance')
ylabel('Pressure, kPa')
legend([h5 h1 h2], 'Empty', 'Sample 1', 'Sample 2')
% legend([h5 h1 h2 h3 h4], 'Empty', 'Sample 1', 'Sample 2', 'Sample 3', 'Sample 4')
set(gca,'FontSize',24)
box on 

% ymin = -10000;
% ymax = 50000;
% figure(Position=[400,500,780,500])
% TL = tiledlayout(2,2,"TileSpacing","tight","Padding",'tight');
% nexttile
% hold on
% plot(normalTime,wdenoise(a25_40p3(sig,:),14),'-b')
% xlim([0 0.1])
% ylim([ymin ymax])
% scatter(normalTime(dataPeakLocs(sig,1:end-2)),dataPeakVals(sig,1:end-2),75,'_r','LineWidth',0.75)
% scatter(normalTime(dataTroughLocs(sig,1:end-1)),dataTroughVals(sig,1:end-1),75,'_r','LineWidth',0.75)
% title('(a)','FontWeight','bold','FontAngle','italic')
% sigAx = gca;
% sigAx.FontSize = 14;
% grid on
% box on 
% 
% nexttile
% hold on
% plot(normalTime,wdenoise(a18_40p3(sig,:),14),'-b')
% xlim([0 0.1])
% ylim([ymin ymax])
% scatter(normalTime(peakLocs2(sig,1:end-2)),peakVals2(sig,1:end-2),75,'_r','LineWidth',0.75)
% scatter(normalTime(troughLocs2(sig,1:end-1)),troughVals2(sig,1:end-1),75,'_r','LineWidth',0.75)
% title('(b)','FontWeight','bold','FontAngle','italic')
% sigAx = gca;
% sigAx.FontSize = 14;
% grid on
% box on 
% 
% nexttile
% hold on
% plot(normalTime,wdenoise(a12_40p3(sig,:),14),'-b')
% xlim([0 0.1])
% ylim([ymin ymax])
% scatter(normalTime(peakLocs3(sig,1:end-2)),peakVals3(sig,1:end-2),75,'_r','LineWidth',0.75)
% scatter(normalTime(troughLocs3(sig,1:end-1)),troughVals3(sig,1:end-1),75,'_r','LineWidth',0.75)
% title('(c)','FontWeight','bold','FontAngle','italic')
% sigAx = gca;
% sigAx.FontSize = 14;
% grid on
% box on 
% 
% nexttile
% hold on
% plot(normalTime,wdenoise(a6_40p3(sig,:),14),'-b')
% xlim([0 0.1])
% ylim([ymin ymax])
% scatter(normalTime(peakLocs4(sig,1:end-2)),peakVals4(sig,1:end-2),75,'_r','LineWidth',0.75)
% scatter(normalTime(troughLocs4(sig,1:end-1)),troughVals4(sig,1:end-1),75,'_r','LineWidth',0.75)
% title('(d)','FontWeight','bold','FontAngle','italic')
% sigAx = gca;
% sigAx.FontSize = 14;
% grid on
% box on 
% 
% xlabel(TL,'Time, s','FontSize',16)
% ylabel(TL,'Pressure, Pa','FontSize',16)
