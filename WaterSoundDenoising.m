clc;
clear all;
close all;

%% Download Audio -------------------------------------------------------------------
disp('Downloading your voice file...');
url = 'https://raw.githubusercontent.com/Himel-Sarder/DSP_Project/main/NoisyVoice.wav';
filename = 'NoisyVoice.wav';
websave(filename, url);
disp('Download complete!');

%% Read Audio -----------------------------------------------------------------------
[y, Fs] = audioread(filename);
t = (0:length(y)-1)/Fs;

%% Play Original Audio --------------------------------------------------------------
disp("Playing Original Voice...");
sound(y, Fs);
pause(length(y)/Fs + 1);

%% Plot Original Waveform -----------------------------------------------------------
figure;
plot(t, y, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('Waveform of Original Voice');
grid on;

%% ---------------------- LOW-PASS FILTER ----------------------
disp('Applying Low-pass Filter...');
Fc_low = 1500;   % Cutoff frequency
[b_low, a_low] = butter(6, Fc_low/(Fs/2), 'low');
y_low = filter(b_low, a_low, y);

disp('Playing Low-pass Filtered Voice...');
sound(y_low, Fs);
pause(length(y)/Fs + 1);

%% ---------------------- HIGH-PASS FILTER ----------------------
disp('Applying High-pass Filter...');
Fc_high = 1000;  % Cutoff frequency
[b_high, a_high] = butter(6, Fc_high/(Fs/2), 'high');
y_high = filter(b_high, a_high, y);

disp('Playing High-pass Filtered Voice...');
sound(y_high, Fs);
pause(length(y)/Fs + 1);

%% ---------------------- BAND-PASS FILTER ----------------------
disp('Applying Band-pass Filter...');
Fc1 = 500; Fc2 = 2000;
[b_band, a_band] = butter(6, [Fc1 Fc2]/(Fs/2), 'bandpass');
y_band = filter(b_band, a_band, y);

disp('Playing Band-pass Filtered Voice...');
sound(y_band, Fs);
pause(length(y)/Fs + 1);

%% ---------------------- DENOISING (LOW-PASS) ----------------------
disp('Applying Noise Reduction...');

% Low-pass to remove high-frequency noise
Fc_denoise = 600;  
[b_denoise, a_denoise] = butter(6, Fc_denoise/(Fs/2), 'low');
y_temp = filter(b_denoise, a_denoise, y);

% Wiener filter (adaptive noise removal)
y_denoised = wiener2(y_temp, [5 1]);

disp('Playing Denoised Voice...');
sound(y_denoised, Fs);
pause(length(y)/Fs + 1);

%% ---------------------- PLOTTING RESULTS ----------------------
figure('Name', 'Audio Filtering and Denoising', 'NumberTitle', 'off');
subplot(5,1,1);
plot(t, y, 'r');
title('Original Audio'); ylabel('Amplitude');
xlim([0 max(t)]);

subplot(5,1,2);
plot(t, y_low, 'color', [1 0.5 0]);
title('Low-pass Filtered'); ylabel('Amplitude');
xlim([0 max(t)]);

subplot(5,1,3);
plot(t, y_high, 'color', [1 0.5 0]);
title('High-pass Filtered'); ylabel('Amplitude');
xlim([0 max(t)]);

subplot(5,1,4);
plot(t, y_band, 'color', [1 0.5 0]);
title('Band-pass Filtered'); ylabel('Amplitude');
xlim([0 max(t)]);

subplot(5,1,5);
plot(t, y_denoised, 'g');
title('Denoised Audio'); xlabel('Time (s)'); ylabel('Amplitude');
xlim([0 max(t)]);

sgtitle('Audio Filtering and Noise Reduction');



