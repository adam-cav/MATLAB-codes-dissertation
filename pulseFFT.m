function [FFTs, normalisedFFTs, NUFFTs, normalisedNUFFTs, freq] = pulseFFT(data,pulseNum,peaksPos,troughsPos,Fs)

L = Fs;
freq = Fs / L * (0:(L / 2));

[FFTs, normalisedFFTs, NUFFTs, normalisedNUFFTs] = deal(zeros(size(data,1), (L/2)+1));
for i = 1:size(data,1)
    lengthCorrection = zeros(1,L);
    pulse = [data(i,peaksPos(i,pulseNum):troughsPos(i,pulseNum)), lengthCorrection(1:end-numel(data(i,peaksPos(i,pulseNum):troughsPos(i,pulseNum))))];
    
    pulseFFT = fft(pulse);

    P2 = abs(pulseFFT / L);
    pulseFFT = P2(1:L / 2 + 1);
    pulseFFT(2:end - 1) = 2 * pulseFFT(2:end - 1);

    FFTs(i,:) = pulseFFT;
    normalisedFFTs(i,:) = pulseFFT/max(abs(pulseFFT));

    t = linspace(0, L / Fs, L);
    f = (0:L-1) / L * Fs;
    pulseNFFT = nufft(pulse, t, f);
    P4 = abs(pulseNFFT / L);
    pulseNFFT = P4(1:L / 2 + 1);
    pulseNFFT(2:end - 1) = 2 * pulseNFFT(2:end - 1);
    
    NUFFTs(i,:) = pulseNFFT;
    normalisedNUFFTs(i,:) = pulseNFFT / max(abs(pulseNFFT));

end