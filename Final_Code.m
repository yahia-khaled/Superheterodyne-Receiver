clc 
clear 
close all
% Load audio file
audioFilePath_1 = 'Short_BBCArabic2.wav';
audioFilePath_2 = 'Short_FM9090.wav';
[audioSignal_BBC, sampleRate_BBC] = audioread(audioFilePath_1);
[audioSignal_FM9090, sampleRate_FM9090] = audioread(audioFilePath_2);
length_BBC=length(audioSignal_BBC);
length_FM9090=length(audioSignal_FM9090);
Orignal_Fs = sampleRate_BBC;
% Extract channels
channel1_BBC = audioSignal_BBC(:, 1);
channel2_BBC = audioSignal_BBC(:, 2);
channel1_FM9090 = audioSignal_FM9090(:, 1);
channel2_FM9090 = audioSignal_FM9090(:, 2);
% Add channels
sumChannels_BBC = channel1_BBC + channel2_BBC;
sumChannels_FM9090 = channel1_FM9090 + channel2_FM9090;
% Plot original signals in Time domain
t_BBC = (0:length(sumChannels_BBC)-1) / sampleRate_BBC; 
t_FM9090 = (0:length(sumChannels_FM9090)-1) / sampleRate_FM9090;  
figure;
plot(t_BBC, sumChannels_BBC);
title('BBC Signal');
xlabel('Time (seconds)');
figure;
plot(t_FM9090, sumChannels_FM9090);
title('FM9090 Signal');
xlabel('Time (seconds)');
% Parameters for FFT
N_BBC = length(sumChannels_BBC(:, 1));
N_FM9090 = length(sumChannels_FM9090(:, 1));
frequencies_BBC = linspace(-sampleRate_BBC/2, sampleRate_BBC/2, N_BBC);
frequencies_FM9090 = linspace(-sampleRate_FM9090/2, sampleRate_FM9090/2, N_FM9090);
fft_BBC = fftshift(fft(sumChannels_BBC(:, 1)));
fft_FM9090= fftshift(fft(sumChannels_FM9090(:, 1)));
% Plot original signals in Frequency domain
figure;
plot(frequencies_BBC, abs(fft_BBC));
title('BBC Signal in Frequency Domain');
xlabel('Frequency (Hz)');
grid on;
figure;
plot(frequencies_FM9090, abs(fft_FM9090));
title('FM9090 Signal in Frequency Domain');
xlabel('Frequency (Hz)');
grid on;
% Modulation
n = 0;          % Signal index
deltaF = 55e3;  % Frequency spacing
carrierFreq = 100e3 + n * deltaF;  % Carrier frequency for the first signal
% Generate carrier signal
if carrierFreq > sampleRate_BBC / 2
    % Increase the sampling frequency
    sampleRate_BBC = sampleRate_BBC * 10;
    sumChannels_BBC = interp(sumChannels_BBC, 10);
    t_BBC = (0:length(sumChannels_BBC)-1) / sampleRate_BBC;
end
carrier = cos(2 * pi * carrierFreq * t_BBC);
% Modulate signals
modulated_BBC = sumChannels_BBC .* carrier';
n = 1;          % Signal index
carrierFreq2 = 100e3 + n * deltaF;  % Carrier frequency for the first signal
% Generate carrier signal
if carrierFreq2 > sampleRate_FM9090 / 2
    % Increase the sampling frequency
    sampleRate_FM9090 = sampleRate_FM9090 * 10;
    sumChannels_FM9090 = interp(sumChannels_FM9090, 10);
    t_FM9090 = (0:length(sumChannels_FM9090)-1) / sampleRate_FM9090;
end
carrier = cos(2 * pi * carrierFreq2 * t_FM9090);
% Modulate signals
modulated_FM9090 = sumChannels_FM9090 .* carrier';
% Parameters for FFT
N_BBC_modulated = length(modulated_BBC(:, 1));
N_FM9090_modulated = length(sumChannels_FM9090(:, 1));
frequencies_BBC = linspace(-sampleRate_BBC/2, sampleRate_BBC/2, N_BBC_modulated);
frequencies_FM9090 = linspace(-sampleRate_FM9090/2, sampleRate_FM9090/2, N_FM9090_modulated);
fft_BBC_modulated = fftshift(fft(modulated_BBC(:, 1)));
fft_FM9090_modulated= fftshift(fft(modulated_FM9090(:, 1)));
% Plot modulated signals in Frequency domain
figure;
plot(frequencies_BBC, abs(fft_BBC_modulated));
title('BBC Modulated Signal in Frequency Domain');
xlabel('Frequency (Hz)');
grid on;
figure;
plot(frequencies_FM9090, abs(fft_FM9090_modulated));
title('FM9090 Modulated Signal in Frequency Domain');
xlabel('Frequency (Hz)');
grid on;
% Ensure both vectors have the same length
maxLength = max(length(modulated_FM9090), length(modulated_BBC));

% Pad the shorter vector with zeros
modulated_FM9090 = [modulated_FM9090; zeros(maxLength - length(modulated_FM9090), 1)];
modulated_BBC = [modulated_BBC; zeros(maxLength - length(modulated_BBC), 1)];

% Now, both vectors have the same length, and you can add them
result_BBC_FM_Summing = modulated_FM9090 + modulated_BBC;
% Parameters for FFT
N_result_BBC_FM_Summing = length(result_BBC_FM_Summing(:, 1));
frequencies_result_BBC_FM_Summing = linspace(-sampleRate_BBC/2, sampleRate_BBC/2, N_result_BBC_FM_Summing);
fft_result_BBC_FM_Summing = fftshift(fft(result_BBC_FM_Summing(:, 1)));
% Plot the Summing Modulated signals in Frequency domain
figure;
plot(frequencies_result_BBC_FM_Summing, abs(fft_result_BBC_FM_Summing));
title('Summing Modulated Signals in Frequency Domain');
xlabel('Frequency (Hz)');
grid on;
%% RF stage
BandPassSpecObj = fdesign.bandpass(91000,93000,107000,110000,80,1,80,441000);
BandPassFilt_100k = design(BandPassSpecObj, 'butter');
fvtool(BandPassFilt_100k) %plot the filter magnitude response
BandPassSpecObj = fdesign.bandpass(14000,148000,161600,170000,80,1,80,441000);
BandPassFilt_150k = design(BandPassSpecObj, 'butter');
fvtool(BandPassFilt_150k) %plot the filter magnitude response
% Apply the 100k Bandpass Filter
BBC_Filtered = filter(BandPassFilt_100k, result_BBC_FM_Summing);
% Apply the 150k Bandpass Filter
FM9090_Filtered = filter(BandPassFilt_150k, result_BBC_FM_Summing);
% Parameters for FFT
N_BBC_Filtered = length(BBC_Filtered(:, 1));
frequencies = linspace(-sampleRate_BBC/2, sampleRate_BBC/2, N_BBC_Filtered);
fft_BBC_Filtered = fftshift(fft(BBC_Filtered(:, 1)));
fft_FM9090_Filtered= fftshift(fft(FM9090_Filtered(:, 1)));
% Plot filtered signals
figure;
plot(frequencies,abs(fft_BBC_Filtered));
title('The BBC signal after using BPF');
xlabel('Frequency (Hz)');
grid on;
figure;
plot(frequencies,abs(fft_FM9090_Filtered));
title('The FM9090 signal after using BPF');
xlabel('Frequency (Hz)');
grid on;
% Show the plots
% Plot the original and filtered signals
figure;
subplot(3, 1, 1);
plot(frequencies_result_BBC_FM_Summing, abs(fft_result_BBC_FM_Summing));
title('Modulated Signal before use filter');
xlabel('Frequency (Hz)');
grid on;

subplot(3, 1, 2);
plot(frequencies,abs(fft_BBC_Filtered));
title('The BBC signal after using BPF');
xlabel('Frequency (Hz)');
grid on;

subplot(3, 1, 3);
plot(frequencies,abs(fft_FM9090_Filtered));
title('The FM9090 signal after using BPF');
xlabel('Frequency (Hz)');
grid on;
%% IF statge
IF_Freq=27.5e3;
offset = 1e3;
IF_Frequency_BBC=carrierFreq + IF_Freq+ offset;
IF_Frequency_FM9090=carrierFreq2 + IF_Freq + offset;
carrier = cos(2 * pi * IF_Frequency_BBC * t_BBC);
IF_BBC = BBC_Filtered .* carrier';
carrier = cos(2 * pi * IF_Frequency_FM9090 * t_BBC);
IF_FM9090 = FM9090_Filtered .* carrier';
% Parameters for FFT
fft_IF_BBC = fftshift(fft(IF_BBC(:, 1)));
fft_IF_FM9090= fftshift(fft(IF_FM9090(:, 1)));
% Plot the demodulation IF
figure;
subplot(2, 1, 1);
plot(frequencies, abs(fft_IF_BBC));
title('BBC signal after demodulation');
xlabel('Frequency (Hz)');
grid on;

subplot(2, 1, 2);
plot(frequencies,abs(fft_IF_FM9090));
title('FM9090 signal after demodulation');
xlabel('Frequency (Hz)');
grid on;
% Using Filtrer
BandPassSpecObj = fdesign.bandpass(18000,20000,35000,38000,80,1,80,441000);
BandPassFilt_IF_BBC = design(BandPassSpecObj, 'butter');
fvtool(BandPassFilt_IF_BBC) %plot the filter magnitude response
BandPassSpecObj = fdesign.bandpass(18500,21000,34000,36500,80,1,80,441000);
BandPassFilt_IF_FM9090 = design(BandPassSpecObj, 'butter');
fvtool(BandPassFilt_IF_FM9090) %plot the filter magnitude response
% Apply The BBF to IF BBC signal
BBC_IF = filter(BandPassFilt_IF_BBC, IF_BBC);
% Apply The BBF to IF FM9090 signal
FM_IF = filter(BandPassFilt_IF_FM9090, IF_FM9090);

fft_BBC_IF= fftshift(fft(BBC_IF(:, 1)));
fft_FM_IF = fftshift(fft(FM_IF(:, 1)));
% Plot the demodulation IF
figure;
subplot(2, 1, 1);
plot(frequencies, abs(fft_BBC_IF));
title('After Use BBF to select BBC signal in IF Freq ');
xlabel('Frequency (Hz)');
grid on;

subplot(2, 1, 2);
plot(frequencies,abs(fft_FM_IF));
title('After Use BBF to select FM9090 signal in IF Freq ');
xlabel('Frequency (Hz)');
grid on;
%% The BaseBand Detection
carrier = cos(2 * pi * IF_Freq * t_BBC);
BaseBand_BBC = BBC_IF .* carrier';
BaseBand_FM9090 = FM_IF .* carrier';
fft_BaseBand_BBC = fftshift(fft(BaseBand_BBC(:, 1)));
fft_BaseBand_FM9090= fftshift(fft(BaseBand_FM9090(:, 1)));
% Plot the demodulation IF
figure;
subplot(2, 1, 1);
plot(frequencies, abs(fft_BaseBand_BBC));
title('BBC signal after Mixing with IF Freq');
xlabel('Frequency (Hz)');
grid on;

subplot(2, 1, 2);
plot(frequencies,abs(fft_BaseBand_FM9090));
title('FM9090 signal after Mixing with IF Freq');
xlabel('Frequency (Hz)');
grid on;
%% Using LowBass filter 
BandPassSpecObj = fdesign.lowpass(7000,8000,1,80,441000);
LowPassFilt_BBC = design(BandPassSpecObj, 'butter');
fvtool(LowPassFilt_BBC) %plot the filter magnitude response
BandPassSpecObj = fdesign.lowpass(6500,7000,1,80,441000);
LowPassFilt_FM9090 = design(BandPassSpecObj, 'butter');
fvtool(LowPassFilt_FM9090) %plot the filter magnitude response
% Apply The BBF to IF BBC signal
BBC_Recieving = filter(LowPassFilt_BBC, BaseBand_BBC);
% Apply The BBF to IF FM9090 signal
FM_Recieving  = filter(LowPassFilt_FM9090, BaseBand_FM9090);

fft_BBC_Recieving = fftshift(fft(BBC_Recieving(:, 1)));
fft_FM_Recieving= fftshift(fft(FM_Recieving(:, 1)));
% Plot the demodulation IF
figure;
subplot(2, 1, 1);
plot(frequencies, abs(fft_BBC_Recieving));
title('BBC signal after Using LowPass filter');
xlabel('Frequency (Hz)');
grid on;

subplot(2, 1, 2);
plot(frequencies,abs(fft_FM_Recieving));
title('FM9090 signal after Using LowPass Filter ');
xlabel('Frequency (Hz)');
grid on;
%% DownSampling
Final_BBC_Audio = downsample(BBC_Recieving,sampleRate_BBC/Orignal_Fs);
Final_FM9090_Audio = downsample(FM_Recieving,sampleRate_BBC/Orignal_Fs);
% Play the modulated signal
sound(Final_BBC_Audio, Orignal_Fs); % Assuming sampleRate_BBC is the sampling frequency
% Wait for the duration of the first audio
pause(length(Final_BBC_Audio) / Orignal_Fs);
sound(Final_FM9090_Audio, Orignal_Fs); % Assuming sampleRate_BBC is the sampling frequency
%% Export Final Signals
audiowrite("BBC_Final_Recieving.wav",Final_BBC_Audio,Orignal_Fs);
audiowrite("FM9090_Final_Recieving.wav",Final_FM9090_Audio,Orignal_Fs);