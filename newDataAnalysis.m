clear 
clc
close all

% e23_0001 and 2: AI1-3 = pos.1-3
% distance between pos 2 and 3 = 1370mm
% distance from pos 3 to end = 178.5
%%%%%%%%%%%%% change troughs to be first 0 cross!! 

%% Constants

gamma  = 1.401; % Specific heat ratio at 20C
P0 = 1.013*10^5; % atmospheric pressure at 20C [Pa]
% p = logspace(log10(1),log10(100000),99999); % Pressure array [Pa]
p = linspace(1,80000,79999);
gasConstR = 287.05; % Specific gas constant of air [J/kg K]
temperature = 20 + 273.15; % Air temperature [k]
rho0 = P0 / (gasConstR * temperature); % density of air
c0 = 343; % sound speed at 20C [m/s]
c = c0 + ((gamma - 1) / 2) * (c0 / P0) .* p; % pressure dependent sound speed [m/s]

% pos.1 
% distanceChange = [6800, 3800];

% pos.2
% distanceChange = [6280, 4320];

% pos.3
distanceChange = [3550, 7050];
%% Signal loading

Fs = 200000;
path = '/Users/adamcavanagh/Desktop/Uni/Year 3/Project/Data/exports_090424/txt/';
transducers = 2:4;

signalLength = round(0.111 * Fs); % Signal length after loading
signalCount = 5; % Number of ruptures 
samples2firstPulse = 750; % number of samples till first pulse

name_e40 = 'e40_000'; % 40 um mylar
name_e23 = 'e23_000'; % 23 um mylar
name_ePa = 'ePa_000'; % single layer paper

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

[e40p1, e40t1, e40p2, e40t2, e40p3, e40t3] = loadMultiAI(name_e40,path,transducers,signalLength,signalCount,samples2firstPulse);
[e23p1, e23t1, e23p2, e23t2, e23p3, e23t3] = loadMultiAI(name_e23,path,transducers,signalLength,signalCount,samples2firstPulse);
[ePap1, ePat1, ePap2, ePat2, ePap3, ePat3] = loadMultiAI(name_ePa,path,transducers,signalLength,signalCount,samples2firstPulse);

[a6_40p1, a6_40t1, a6_40p2, a6_40t2, a6_40p3, a6_40t3] = loadMultiAI(name_a6_40,path,transducers,signalLength,signalCount,samples2firstPulse);
[a6_23p1, a6_23t1, a6_23p2, a6_23t2, a6_23p3, a6_23t3] = loadMultiAI(name_a6_23,path,transducers,signalLength,signalCount,samples2firstPulse);
[a6_Pap1, a6_Pat1, a6_Pap2, a6_Pat2, a6_Pap3, a6_Pat3] = loadMultiAI(name_a6_Pa,path,transducers,signalLength,signalCount,samples2firstPulse);

[a12_40p1, a12_40t1, a12_40p2, a12_40t2, a12_40p3, a12_40t3] = loadMultiAI(name_a12_40,path,transducers,signalLength,signalCount,samples2firstPulse);
[a12_23p1, a12_23t1, a12_23p2, a12_23t2, a12_23p3, a12_23t3] = loadMultiAI(name_a12_23,path,transducers,signalLength,signalCount,samples2firstPulse);
[a12_Pap1, a12_Pat1, a12_Pap2, a12_Pat2, a12_Pap3, a12_Pat3] = loadMultiAI(name_a12_Pa,path,transducers,signalLength,signalCount,samples2firstPulse);

[a18_40p1, a18_40t1, a18_40p2, a18_40t2, a18_40p3, a18_40t3] = loadMultiAI(name_a18_40,path,transducers,signalLength,signalCount,samples2firstPulse);
[a18_23p1, a18_23t1, a18_23p2, a18_23t2, a18_23p3, a18_23t3] = loadMultiAI(name_a18_23,path,transducers,signalLength,signalCount,samples2firstPulse);
[a18_Pap1, a18_Pat1, a18_Pap2, a18_Pat2, a18_Pap3, a18_Pat3] = loadMultiAI(name_a18_Pa,path,transducers,signalLength,signalCount,samples2firstPulse);

[a25_40p1, a25_40t1, a25_40p2, a25_40t2, a25_40p3, a25_40t3] = loadMultiAI(name_a25_40,path,transducers,signalLength,signalCount,samples2firstPulse);
[a25_23p1, a25_23t1, a25_23p2, a25_23t2, a25_23p3, a25_23t3] = loadMultiAI(name_a25_23,path,transducers,signalLength,signalCount,samples2firstPulse);
[a25_Pap1, a25_Pat1, a25_Pap2, a25_Pat2, a25_Pap3, a25_Pat3] = loadMultiAI(name_a25_Pa,path,transducers,signalLength,signalCount,samples2firstPulse);


%% Analysis
% sig = 1;
micPos = 3;

numPeaks = 5; % Number of peaks to extract
minPeakDistance = 2000; % Minimum distance between pulse instances
peakErrorCheck = 500; % Number of samples to error check either side of initial values

%% Empty Tube
micPos = 1;

dataE40 = sprintf('e40p%d',micPos);
timeE40 = sprintf('e40t%d',micPos);
dataE23 = sprintf('e23p%d',micPos);
timeE23 = sprintf('e23t%d',micPos);
dataEPa = sprintf('ePap%d',micPos);
timeEPa = sprintf('ePat%d',micPos);

[e40peaksVals, e40peaksPos, e40troughsVals, e40troughsPos] = extractPeaks(eval(dataE40),eval(timeE40),numPeaks,minPeakDistance,peakErrorCheck,distanceChange);
[e23peaksVals, e23peaksPos, e23troughsVals, e23troughsPos] = extractPeaks(eval(dataE23),eval(timeE23),numPeaks,minPeakDistance,peakErrorCheck,distanceChange);
[ePapeaksVals, ePapeaksPos, ePatroughsVals, ePatroughsPos] = extractPeaks(eval(dataEPa),eval(timeEPa),numPeaks,minPeakDistance,peakErrorCheck,distanceChange);

[e40durationModel, e40pressureModel, distanceVector, e40measuredDurations, measuredDistanceVector, e40modelPressures] = pulseDistanceRelations(eval(dataE40),eval(timeE40),distanceChange,e40peaksVals,e40peaksPos,e40troughsPos);
[e23durationModel, e23pressureModel, ~, e23measuredDurations, ~, e23modelPressures] = pulseDistanceRelations(eval(dataE23),eval(timeE23),distanceChange,e23peaksVals,e23peaksPos,e23troughsPos);
[ePadurationModel, ePapressureModel, ~, ePameasuredDurations, ~, ePamodelPressures] = pulseDistanceRelations(eval(dataEPa),eval(timeEPa),distanceChange,ePapeaksVals,ePapeaksPos,ePatroughsPos);

%% ABH 6 cavities
% micPos = 1;

data640 = sprintf('a6_40p%d',micPos);
time640 = sprintf('a6_40t%d',micPos);
data623 = sprintf('a6_23p%d',micPos);
time623 = sprintf('a6_23t%d',micPos);
data6Pa = sprintf('a6_Pap%d',micPos);
time6Pa = sprintf('a6_Pat%d',micPos);

[a6_40peaksVals, a6_40peaksPos, a6_40troughsVals, a6_40troughsPos] = extractPeaks(eval(data640),eval(time640),numPeaks,minPeakDistance,peakErrorCheck,distanceChange);
[a6_23peaksVals, a6_23peaksPos, a6_23troughsVals, a6_23troughsPos] = extractPeaks(eval(data623),eval(time623),numPeaks,minPeakDistance,peakErrorCheck,distanceChange);
[a6_PapeaksVals, a6_PapeaksPos, a6_PatroughsVals, a6_PatroughsPos] = extractPeaks(eval(data6Pa),eval(time6Pa),numPeaks,minPeakDistance,peakErrorCheck,distanceChange);

[a6_40durationModel, a6_40pressureModel, ~, a6_40measuredDurations, ~, a6_40modelPressures] = pulseDistanceRelations(eval(data640),eval(time640),distanceChange,a6_40peaksVals,a6_40peaksPos,a6_40troughsPos);
[a6_23durationModel, a6_23pressureModel, ~, a6_23measuredDurations, ~, a6_23modelPressures] = pulseDistanceRelations(eval(data623),eval(time623),distanceChange,a6_23peaksVals,a6_23peaksPos,a6_23troughsPos);
[a6_PadurationModel, a6_PapressureModel, ~, a6_PameasuredDurations, ~, a6_PamodelPressures] = pulseDistanceRelations(eval(data6Pa),eval(time6Pa),distanceChange,a6_PapeaksVals,a6_PapeaksPos,a6_PatroughsPos);

%% ABH 12 cavities
% micPos = 1;

data1240 = sprintf('a12_40p%d',micPos);
time1240 = sprintf('a12_40t%d',micPos);
data1223 = sprintf('a12_23p%d',micPos);
time1223 = sprintf('a12_23t%d',micPos);
data12Pa = sprintf('a12_Pap%d',micPos);
time12Pa = sprintf('a12_Pat%d',micPos);

[a12_40peaksVals, a12_40peaksPos, a12_40troughsVals, a12_40troughsPos] = extractPeaks(eval(data1240),eval(time1240),numPeaks,minPeakDistance,peakErrorCheck,distanceChange);
[a12_23peaksVals, a12_23peaksPos, a12_23troughsVals, a12_23troughsPos] = extractPeaks(eval(data1223),eval(time1223),numPeaks,minPeakDistance,peakErrorCheck,distanceChange);
[a12_PapeaksVals, a12_PapeaksPos, a12_PatroughsVals, a12_PatroughsPos] = extractPeaks(eval(data12Pa),eval(time12Pa),numPeaks,minPeakDistance,peakErrorCheck,distanceChange);

[a12_40durationModel, a12_40pressureModel, ~, a12_40measuredDurations, ~, a12_40modelPressures] = pulseDistanceRelations(eval(data1240),eval(time1240),distanceChange,a12_40peaksVals,a12_40peaksPos,a12_40troughsPos);
[a12_23durationModel, a12_23pressureModel, ~, a12_23measuredDurations, ~, a12_23modelPressures] = pulseDistanceRelations(eval(data1223),eval(time1223),distanceChange,a12_23peaksVals,a12_23peaksPos,a12_23troughsPos);
[a12_PadurationModel, a12_PapressureModel, ~, a12_PameasuredDurations, ~, a12_PamodelPressures] = pulseDistanceRelations(eval(data12Pa),eval(time12Pa),distanceChange,a12_PapeaksVals,a12_PapeaksPos,a12_PatroughsPos);

%% ABH 18 cavities
% micPos = 1;

data1840 = sprintf('a18_40p%d',micPos);
time1840 = sprintf('a18_40t%d',micPos);
data1823 = sprintf('a18_23p%d',micPos);
time1823 = sprintf('a18_23t%d',micPos);
data18Pa = sprintf('a18_Pap%d',micPos);
time18Pa = sprintf('a18_Pat%d',micPos);

[a18_40peaksVals, a18_40peaksPos, a18_40troughsVals, a18_40troughsPos] = extractPeaks(eval(data1840),eval(time1840),numPeaks,minPeakDistance,peakErrorCheck,distanceChange);
[a18_23peaksVals, a18_23peaksPos, a18_23troughsVals, a18_23troughsPos] = extractPeaks(eval(data1823),eval(time1823),numPeaks,minPeakDistance,peakErrorCheck,distanceChange);
[a18_PapeaksVals, a18_PapeaksPos, a18_PatroughsVals, a18_PatroughsPos] = extractPeaks(eval(data18Pa),eval(time18Pa),numPeaks,minPeakDistance,peakErrorCheck,distanceChange);

[a18_40durationModel, a18_40pressureModel, ~, a18_40measuredDurations, ~, a18_40modelPressures] = pulseDistanceRelations(eval(data1840),eval(time1840),distanceChange,a18_40peaksVals,a18_40peaksPos,a18_40troughsPos);
[a18_23durationModel, a18_23pressureModel, ~, a18_23measuredDurations, ~, a18_23modelPressures] = pulseDistanceRelations(eval(data1823),eval(time1823),distanceChange,a18_23peaksVals,a18_23peaksPos,a18_23troughsPos);
[a18_PadurationModel, a18_PapressureModel, ~, a18_PameasuredDurations, ~, a18_PamodelPressures] = pulseDistanceRelations(eval(data18Pa),eval(time18Pa),distanceChange,a18_PapeaksVals,a18_PapeaksPos,a18_PatroughsPos);

%% ABH 25 cavities
% micPos = 1;

data2540 = sprintf('a25_40p%d',micPos);
time2540 = sprintf('a25_40t%d',micPos);
data2523 = sprintf('a25_23p%d',micPos);
time2523 = sprintf('a25_23t%d',micPos);
data25Pa = sprintf('a25_Pap%d',micPos);
time25Pa = sprintf('a25_Pat%d',micPos);

[a25_40peaksVals, a25_40peaksPos, a25_40troughsVals, a25_40troughsPos] = extractPeaks(eval(data2540),eval(time2540),numPeaks,minPeakDistance,peakErrorCheck,distanceChange);
[a25_23peaksVals, a25_23peaksPos, a25_23troughsVals, a25_23troughsPos] = extractPeaks(eval(data2523),eval(time2523),numPeaks,minPeakDistance,peakErrorCheck,distanceChange);
[a25_PapeaksVals, a25_PapeaksPos, a25_PatroughsVals, a25_PatroughsPos] = extractPeaks(eval(data25Pa),eval(time25Pa),numPeaks,minPeakDistance,peakErrorCheck,distanceChange);

[a25_40durationModel, a25_40pressureModel, ~, a25_40measuredDurations, ~, a25_40modelPressures] = pulseDistanceRelations(eval(data2540),eval(time2540),distanceChange,a25_40peaksVals,a25_40peaksPos,a25_40troughsPos);
[a25_23durationModel, a25_23pressureModel, ~, a25_23measuredDurations, ~, a25_23modelPressures] = pulseDistanceRelations(eval(data2523),eval(time2523),distanceChange,a25_40peaksVals,a25_40peaksPos,a25_40troughsPos);
[a25_PadurationModel, a25_PapressureModel, ~, a25_PameasuredDurations, ~, a25_PamodelPressures] = pulseDistanceRelations(eval(data25Pa),eval(time25Pa),distanceChange,a25_PapeaksVals,a25_PapeaksPos,a25_PatroughsPos);

%% Energies

a25_40pulseEnergyMatrix = pulseEnergies(eval(data2540),a25_40peaksPos,a25_40troughsPos);
a18_40pulseEnergyMatrix = pulseEnergies(eval(data1840),a18_40peaksPos,a18_40troughsPos);
a12_40pulseEnergyMatrix = pulseEnergies(eval(data1240),a12_40peaksPos,a12_40troughsPos);
a6_40pulseEnergyMatrix = pulseEnergies(eval(data640),a6_40peaksPos,a6_40troughsPos);

a25_23pulseEnergyMatrix = pulseEnergies(eval(data2523),a25_23peaksPos,a25_23troughsPos);
a18_23pulseEnergyMatrix = pulseEnergies(eval(data1823),a18_23peaksPos,a18_23troughsPos);
a12_23pulseEnergyMatrix = pulseEnergies(eval(data1223),a12_23peaksPos,a12_23troughsPos);
a6_23pulseEnergyMatrix = pulseEnergies(eval(data623),a6_23peaksPos,a6_23troughsPos);

a25_PapulseEnergyMatrix = pulseEnergies(eval(data25Pa),a25_PapeaksPos,a25_PatroughsPos);
a18_PapulseEnergyMatrix = pulseEnergies(eval(data18Pa),a18_PapeaksPos,a18_PatroughsPos);
a12_PapulseEnergyMatrix = pulseEnergies(eval(data12Pa),a12_PapeaksPos,a12_PatroughsPos);
a6_PapulseEnergyMatrix = pulseEnergies(eval(data6Pa),a6_PapeaksPos,a6_PatroughsPos);

e40_pulseEnergyMatrix = pulseEnergies(eval(dataE40),e40peaksPos,e40troughsPos);
e23_pulseEnergyMatrix = pulseEnergies(eval(dataE23),e23peaksPos,e23troughsPos);
ePa_pulseEnergyMatrix = pulseEnergies(eval(dataEPa),ePapeaksPos,ePatroughsPos);

[~,a25energyAbsCoeff] = energyReflAbsCoefficient(eval(data2540),a25_40pulseEnergyMatrix,a25_40peaksPos);
[~,a18energyAbsCoeff] = energyReflAbsCoefficient(eval(data1840),a18_40pulseEnergyMatrix,a18_40peaksPos);
[~,a12energyAbsCoeff] = energyReflAbsCoefficient(eval(data1240),a12_40pulseEnergyMatrix,a12_40peaksPos);
[~,a6energyAbsCoeff] = energyReflAbsCoefficient(eval(data640),a6_40pulseEnergyMatrix,a6_40peaksPos);

%% Averaging

e40durationModel = mean(e40durationModel,1);                e40pressureModel = mean(e40pressureModel,1);
e23durationModel = mean(e23durationModel,1);                e23pressureModel = mean(e23pressureModel,1);
ePadurationModel = mean(ePadurationModel,1);                ePapressureModel = mean(ePapressureModel,1);
e40measuredDurations = mean(e40measuredDurations,1);        e40modelPressures = mean(e40modelPressures,1);
e23measuredDurations = mean(e23measuredDurations,1);        e23modelPressures = mean(e23modelPressures,1);
ePameasuredDurations = mean(ePameasuredDurations,1);        ePamodelPressures = mean(ePamodelPressures,1);

a6_40durationModel = mean(a6_40durationModel,1);            a6_40pressureModel = mean(a6_40pressureModel,1);
a6_23durationModel = mean(a6_23durationModel,1);            a6_23pressureModel = mean(a6_23pressureModel,1);
a6_PadurationModel = mean(a6_PadurationModel,1);            a6_PapressureModel = mean(a6_PapressureModel,1);
a6_40measuredDurations = mean(a6_40measuredDurations,1);    a6_40modelPressures = mean(a6_40modelPressures,1);
a6_23measuredDurations = mean(a6_23measuredDurations,1);    a6_23modelPressures = mean(a6_23modelPressures,1);
a6_PameasuredDurations = mean(a6_PameasuredDurations,1);    a6_PamodelPressures = mean(a6_PamodelPressures,1);

a12_40durationModel = mean(a12_40durationModel,1);          a12_40pressureModel = mean(a12_40pressureModel,1);
a12_23durationModel = mean(a12_23durationModel,1);          a12_23pressureModel = mean(a12_23pressureModel,1);
a12_PadurationModel = mean(a12_PadurationModel,1);          a12_PapressureModel = mean(a12_PapressureModel,1);
a12_40measuredDurations = mean(a12_40measuredDurations,1);  a12_40modelPressures = mean(a12_40modelPressures,1);
a12_23measuredDurations = mean(a12_23measuredDurations,1);  a12_23modelPressures = mean(a12_23modelPressures,1);
a12_PameasuredDurations = mean(a12_PameasuredDurations,1);  a12_PamodelPressures = mean(a12_PamodelPressures,1);

a18_40durationModel = mean(a18_40durationModel,1);          a18_40pressureModel = mean(a18_40pressureModel,1);
a18_23durationModel = mean(a18_23durationModel,1);          a18_23pressureModel = mean(a18_23pressureModel,1);
a18_PadurationModel = mean(a18_PadurationModel,1);          a18_PapressureModel = mean(a18_PapressureModel,1);
a18_40measuredDurations = mean(a18_40measuredDurations,1);  a18_40modelPressures = mean(a18_40modelPressures,1);
a18_23measuredDurations = mean(a18_23measuredDurations,1);  a18_23modelPressures = mean(a18_23modelPressures,1);
a18_PameasuredDurations = mean(a18_PameasuredDurations,1);  a18_PamodelPressures = mean(a18_PamodelPressures,1);

a25_40durationModel = mean(a25_40durationModel,1);          a25_40pressureModel = mean(a25_40pressureModel,1);
a12_23durationModel = mean(a12_23durationModel,1);          a12_23pressureModel = mean(a12_23pressureModel,1);
a25_PadurationModel = mean(a25_PadurationModel,1);          a25_PapressureModel = mean(a25_PapressureModel,1);
a25_40measuredDurations = mean(a25_40measuredDurations,1);  a25_40modelPressures = mean(a25_40modelPressures,1);
a12_23measuredDurations = mean(a12_23measuredDurations,1);  a12_23modelPressures = mean(a12_23modelPressures,1);
a25_PameasuredDurations = mean(a25_PameasuredDurations,1);  a25_PamodelPressures = mean(a25_PamodelPressures,1);

%% Mach numbersss



%% Plots
close all

% % %% Empty tube distance relationship plot

xtikMax = ceil(measuredDistanceVector(1:numPeaks)/10)*10;
majorTicks = 0:2.5:max(xtikMax);

edistanceFig = figure;
edistanceFig.Position = [902,3,897,1042];
eTL = tiledlayout(2,1,"TileSpacing",'tight','Padding','none');
nexttile 
hold on 
plot(distanceVector,e40durationModel,'--r')
scatter(measuredDistanceVector,e40measuredDurations,40,'ro','filled')
xticks(majorTicks)
xlim([0 35])
grid on 
box on
ylabel('Pulse duration, s','FontSize',18)
eAx1 = gca;

nexttile
hold on
plot(distanceVector,e40pressureModel,'--b')
scatter(measuredDistanceVector,mean(e40peaksVals,1),'bo','filled')
xticks(majorTicks)
xlim([0 35])
grid on 
box on
xlabel(eTL,'Distance travelled, m','FontSize',18)
ylabel('Pressure, Pa','FontSize',18)
eAx2 = gca;
set([eAx1 eAx2], 'FontSize',18)

% %% abh rigid termination plot

sig = 1;
pulse = 1;
offset = 3; % for pulse plots, in samples 
length = 1500;

a25rigidData = a25_23p3;
a25peaks = a25_23peaksPos;
a25peaksVals = a25_23peaksVals;
a25troughs = a25_23troughsPos;
a25troughsVals = a25_23troughsVals;

a18rigidData = a18_Pap3;
a18peaks = a18_PapeaksPos;
a18peaksVals = a18_PapeaksVals;
a18troughs = a18_PatroughsPos;
a18troughsVals = a18_PatroughsVals;

a12rigidData = a12_40p3;
a12peaks = a12_40peaksPos;
a12peaksVals = a12_40peaksVals;
a12troughs = a12_40troughsPos;
a12troughsVals = a12_40troughsVals;

a6rigidData = a6_40p3;
a6peaks = a6_40peaksPos;
a6peaksVals = a6_PapeaksVals;
a6troughs = a6_PatroughsPos;
a6troughsVals = a6_PatroughsVals;

e40data = e40p3;
e40peaks = e40peaksPos;
e40troughs = ePatroughsPos;

figure('Position',[2, 500, 800, 650])
tiledlayout(1,1,"TileSpacing","none","Padding",'none')
nexttile 
hold on 
plot(a6rigidData(sig,a6peaks(sig,pulse)-offset:end),'-k','LineWidth',1)
plot(a6rigidData(sig,a6peaks(sig,pulse+1)-offset:end),'-r','LineWidth',1)
plot(e40data(sig,e40peaks(sig,pulse)-offset:end),'-b','LineWidth',0.75)
plot(e40data(sig,e40peaks(sig,pulse+1)-offset:end),'-','Color',[0.5 0 0.8],'LineWidth',0.75)
xlabel('Time (ms)')
ylabel('Pressure (kPa)')
xticks([0, 100, 200, 300, 400, 500, 600, 700, 800])
xticklabels([0, 100, 200, 300, 400, 500, 600, 700, 800] .* 1000 ./ Fs)
% xticklabels(xticks./Fs)
yticklabels(yticks ./ 1000)
xlim([0 800])
box on 
legend('Incident, ABH','Reflected, ABH','Incident, empty','Reflected, empty')
set(gca,'FontSize',18)

% incReflFig = figure;
% % incReflFig.Position = [3,565,880,700];
% incReflFig.Position = [3,565,897,1024];
% incReflTL = tiledlayout(2,1,"TileSpacing",'tight','Padding','none');
% nexttile 
% hold on 
% plot(a25rigidData(sig,a25peaks(sig,pulse)-offset:a25troughs(sig,pulse)+offset),'-k','LineWidth',0.75)
% plot(e40data(sig,e40peaks(sig,pulse)-offset:e40troughs(sig,pulse)+offset),'-b','LineWidth',0.75)
% plot(e40data(sig,e40peaks(sig,pulse+1)-offset:e40troughs(sig,pulse+1)+offset),'--b','LineWidth',0.75)
% plot(a25rigidData(sig,a25peaks(sig,pulse+1)-offset:a25troughs(sig,pulse+1)+offset),'--k','LineWidth',0.75)
% xlim([0 1000])
% grid on 
% box on
% 
% nexttile
% hold on
% plot(a18rigidData(sig,a18peaks(sig,pulse)-offset:a18troughs(sig,pulse)+offset),'-k','LineWidth',0.75)
% plot(e40data(sig,e40peaks(sig,pulse)-offset:e40troughs(sig,pulse)+offset),'-b','LineWidth',0.75)
% plot(e40data(sig,e40peaks(sig,pulse+1)-offset:e40troughs(sig,pulse+1)+offset),'--b','LineWidth',0.75)
% plot(a18rigidData(sig,a18peaks(sig,pulse+1)-offset:a18troughs(sig,pulse+1)+offset),'-r','LineWidth',0.75)
% xlim([0 1500])
% 
% xlabel(incReflTL,'Samples')
% ylabel(incReflTL,'Pressure, Pa')
% 
% 
% incReflFig1 = figure;
% % incReflFig1.Position = [3,565,880,700];
% incReflFig1.Position = [904,23,897,1024];
% incReflTL1 = tiledlayout(2,1,"TileSpacing",'tight','Padding','none');
% nexttile 
% hold on 
% plot(a12rigidData(sig,a12peaks(sig,pulse)-offset:a12troughs(sig,pulse)+offset),'-k','LineWidth',0.75)
% plot(e40data(sig,e40peaks(sig,pulse)-offset:e40troughs(sig,pulse)+offset),'-b','LineWidth',0.75)
% plot(e40data(sig,e40peaks(sig,pulse+1)-offset:e40troughs(sig,pulse+1)+offset),'--b','LineWidth',0.75)
% plot(a12rigidData(sig,a12peaks(sig,pulse+1)-offset:a12troughs(sig,pulse+1)+offset),'--k','LineWidth',0.75)
% xlim([0 1000])
% grid on 
% box on
% 
% nexttile
% hold on
% plot(a6rigidData(sig,a6peaks(sig,pulse)-offset:a6troughs(sig,pulse)+offset),'-k','LineWidth',0.75)
% plot(e40data(sig,e40peaks(sig,pulse)-offset:e40troughs(sig,pulse)+offset),'-b','LineWidth',0.75)
% plot(e40data(sig,e40peaks(sig,pulse+1)-offset:e40troughs(sig,pulse+1)+offset),'--b','LineWidth',0.75)
% plot(a6rigidData(sig,a6peaks(sig,pulse+1)-offset:a6troughs(sig,pulse+1)+offset),'-r','LineWidth',0.75)
% xlim([0 1500])
% 
% xlabel(incReflTL1,'Samples')
% ylabel(incReflTL1,'Pressure, Pa')

% % %% Peak checks
%%
sig = 1;
data = ePap1;
% data = wdenoise(data,1);
peaks = ePapeaksVals;
peaksPos = ePapeaksPos;
troughs = ePatroughsVals;
troughsPos = ePatroughsPos;

peakCheckFig = figure;
% peakCheckFig.Position = [1,565,750,600];
peakCheckTL = tiledlayout(5,1,"TileSpacing","tight","Padding","none");

colours = {'xk','xr','xb','xg','xy'};

for i = 1:5
    nexttile
    plot(data(i,:),'-k',"LineWidth",0.1)
    hold on
    for j = 1:5
        scatter(peaksPos(i,j),peaks(i,j),string(colours(j)),'LineWidth',2)
        scatter(troughsPos(i,j),troughs(i,j),string(colours(j)),'LineWidth',2)
    end
end


% 
% figure('Position',[620,528,560,440])
% yyaxis right
% semilogx(mean(abs(fft),1),'-r','LineWidth',2.5)
% xlabel('Frequency, f (Hz)',"Interpreter","tex")
% ylabel('Normalised Magnitude')
% ax = gca;
% % xticks([125, 250, 500, 1000, 2000, 4000])
% xlim([100 5000])
% ylim([0 1.4])
% % grid on
% % box off
% % xline([100, 160, 200, 315, 400, 630, 790, 1260, 1580, 2500, 3200, 5000], '-k','Alpha',1,LineWidth=1)
% % ax.XMinorGrid = 'off';
% ax.TickDir = 'none';
% % ax.FontSize = 13;
% % ax.GridAlpha = 1;
% % ax.GridLineWidth = 1;

%% 

[pks, locs, w,p] = findpeaks(data(1,:),Fs);