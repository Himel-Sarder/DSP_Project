clc;
clear all;
close all;

%% ---------------------- DOWNLOAD AUDIO ----------------------
disp('Downloading your voice file...');
url = 'https://raw.githubusercontent.com/Himel-Sarder/DSP_Project/main/WindNoise.wav';
filename = 'WindNoise.wav';
websave(filename, url);
disp('Download complete!');

%% ---------------------- READ AUDIO ----------------------
[y, Fs] = audioread(filename);
t = (0:length(y)-1)/Fs;

%% ---------------------- PLAY ORIGINAL AUDIO ----------------------
disp("Playing Original Audio...");
sound(y, Fs);
pause(length(y)/Fs + 1);

%% ---------------------- PLOT ORIGINAL WAVEFORM ----------------------
figure('Name', 'Original Audio', 'NumberTitle', 'off');
plot(t, y, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('Waveform of Original Audio');
grid on;

%% ---------------------- LOW-PASS FILTER ----------------------
disp('Applying Low-pass Filter...');
Fc_low = 1500;   % Cutoff frequency in Hz
order_low = 6;   % Filter order
[b_low, a_low] = butter(order_low, Fc_low/(Fs/2), 'low');
y_low = filter(b_low, a_low, y);

disp('Playing Low-pass Filtered Audio...');
sound(y_low, Fs);
pause(length(y)/Fs + 1);

%% ---------------------- HIGH-PASS FILTER ----------------------
disp('Applying High-pass Filter...');
Fc_high = 1000;  % Cutoff frequency in Hz
order_high = 6;  % Filter order
[b_high, a_high] = butter(order_high, Fc_high/(Fs/2), 'high');
y_high = filter(b_high, a_high, y);

disp('Playing High-pass Filtered Audio...');
sound(y_high, Fs);
pause(length(y)/Fs + 1);

%% ---------------------- BAND-PASS FILTER ----------------------
disp('Applying Band-pass Filter...');
Fc1 = 500; Fc2 = 2000;  % Frequency range in Hz
order_band = 6;
[b_band, a_band] = butter(order_band, [Fc1 Fc2]/(Fs/2), 'bandpass');
y_band = filter(b_band, a_band, y);

disp('Playing Band-pass Filtered Audio...');
sound(y_band, Fs);
pause(length(y)/Fs + 1);

%% ---------------------- HIGH-PASS FILTER FOR DENOISING ----------------------
disp('Applying High-pass Filter for Denoising...');

Fc_denoise = 1200;  % Cutoff frequency for denoising in Hz
order_denoise = 6;  % Filter order
[b_denoise, a_denoise] = butter(order_denoise, Fc_denoise/(Fs/2), 'high');

% Apply high-pass filter
y_temp = filter(b_denoise, a_denoise, y);

% Wiener filter for adaptive noise removal
y_denoised = wiener2(y_temp, [5 1]);

disp('Playing Denoised Audio...');
sound(y_denoised, Fs);
pause(length(y)/Fs + 1);

%% ---------------------- PLOTTING ALL RESULTS ----------------------
figure('Name', 'Audio Filtering and Denoising', 'NumberTitle', 'off');

subplot(5,1,1);
plot(t, y, 'r'); 
title('Original Audio'); 
ylabel('Amplitude'); 
xlim([0 max(t)]);
grid on;

subplot(5,1,2);
plot(t, y_low, 'color', [1 0.5 0]); 
title('Low-pass Filtered'); 
ylabel('Amplitude'); 
xlim([0 max(t)]);
grid on;

subplot(5,1,3);
plot(t, y_high, 'color', [0 0.5 1]); 
title('High-pass Filtered'); 
ylabel('Amplitude'); 
xlim([0 max(t)]);
grid on;

subplot(5,1,4);
plot(t, y_band, 'color', [0.5 0 0.5]); 
title('Band-pass Filtered'); 
ylabel('Amplitude'); 
xlim([0 max(t)]);
grid on;

subplot(5,1,5);
plot(t, y_denoised, 'g'); 
title('Denoised Audio'); 
xlabel('Time (s)'); 
ylabel('Amplitude'); 
xlim([0 max(t)]);
grid on;

sgtitle('Audio Filtering and Noise Reduction');
