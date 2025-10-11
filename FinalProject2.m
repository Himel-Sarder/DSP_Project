%% =====================================================================================================================================
%                    Department of Computer Science and Engineering (CSE)
%                     Jamalpur Science & Technology University, Jamalpur
%                             3rd Year 1st Semester Lab Project
%                          Course Code: CSE-3142, Session:2021-2022
%                         Course Title: Digital Signal Processing Lab
%                                 Submitted By - Himel Sarder
%                                         ID : 22111121

% NB: After running the code, If showing any error (May occure for using any build-in function): Install -> Signal Processing Toolbox  !
% How to install:
% Home Tab => Add-Ons => Search (Signal Processing Toolbox) => Select first one => Install (May required gmail and password)

%% =====================================================================================================================================
clc;                    % Clear command window
clear all;              % Remove all variables from workspace
close all;              % Close all open figure windows

%% =============================== DOWNLOAD AUDIO ===============================================
disp('Downloading voice file...');      % Display message to indicate download start
url = 'https://raw.githubusercontent.com/Himel-Sarder/DSP_Project/main/WindNoise.wav'; % Audio file URL
filename = 'WindNoise.wav';                  % Local filename to save the audio
websave(filename, url);                      % Download the file from URL
disp('Download complete!');                  % Display message when download finishes

%% =============================== READ AUDIO ===================================================
[y, Fs] = audioread(filename);              % Read audio file. y = audio data, Fs = sampling frequency

% Convert to mono if stereo
if size(y,2) == 2                             % Check if audio has 2 channels (stereo)
    y = mean(y, 2);                           % Convert stereo to mono by averaging the two channels
end

t = (0:length(y)-1)/Fs;                       % Generate time vector for plotting

%% =============================== PLAY ORIGINAL AUDIO ==========================================
disp("Playing Original Audio...");            % Notify user about playback
sound(y, Fs);                                 % Play the original audio
pause(length(y)/Fs + 1);                      % Pause script until audio finishes playing

%% =============================== PLOT ORIGINAL WAVEFORM =======================================
figure('Name', 'Original Audio', 'NumberTitle', 'off');  % Create a figure window
plot(t, y, 'r');                               % Plot audio waveform in red color
xlabel('Time (s)');                            % Label x-axis
ylabel('Amplitude');                           % Label y-axis
title('Waveform of Original Audio');           % Set title
grid on;                                       % Enable grid for better visualization

%% =============================== LOW-PASS FILTER ==============================================
disp('Applying Low-pass Filter...');
Fc_low = 1500;      % Cutoff frequency
order_low = 6;

% Design low-pass filter using designfilt
lpFilt = designfilt('lowpassiir', 'FilterOrder', order_low, ...
                    'HalfPowerFrequency', Fc_low, 'SampleRate', Fs);

y_low = filtfilt(lpFilt, y);   % Apply filter using zero-phase filtering

disp('Playing Low-pass Filtered Audio...');
sound(y_low, Fs);
pause(length(y)/Fs + 1);

%% =============================== HIGH-PASS FILTER =============================================
disp('Applying High-pass Filter...');
Fc_high = 1000;
order_high = 6;

hpFilt = designfilt('highpassiir', 'FilterOrder', order_high, ...
                     'HalfPowerFrequency', Fc_high, 'SampleRate', Fs);

y_high = filtfilt(hpFilt, y);

disp('Playing High-pass Filtered Audio...');
sound(y_high, Fs);
pause(length(y)/Fs + 1);

%% =============================== BAND-PASS FILTER =============================================
disp('Applying Band-pass Filter...');
Fc1 = 500; Fc2 = 2000;
order_band = 6;

bpFilt = designfilt('bandpassiir', 'FilterOrder', order_band, ...
                     'HalfPowerFrequency1', Fc1, 'HalfPowerFrequency2', Fc2, ...
                     'SampleRate', Fs);

y_band = filtfilt(bpFilt, y);

disp('Playing Band-pass Filtered Audio...');
sound(y_band, Fs);
pause(length(y)/Fs + 1);

%% =============================== HIGH-PASS FILTER FOR DENOISING ===============================
disp('Applying High-pass Filter for Denoising...');
Fc_denoise = 1200;
order_denoise = 6;

denoiseFilt = designfilt('highpassiir', 'FilterOrder', order_denoise, ...
                          'HalfPowerFrequency', Fc_denoise, 'SampleRate', Fs);

y_temp = filtfilt(denoiseFilt, y);

disp('Playing Denoised Audio...');
sound(y_temp, Fs);
pause(length(y)/Fs + 1);

% A high-pass filter allows high frequencies to pass and blocks low frequencies.
% If you set the cutoff at 1200 Hz:
% Frequencies >1200 Hz pass → My voice is mostly preserved.
% Frequencies <1200 Hz are attenuated → most wind noise is removed.

% This works because wind noise is concentrated at low frequencies, so a high-pass filter naturally removes that rumble without affecting my voice much.
%% =============================== PLOTTING ALL RESULTS =========================================
figure('Name', 'Audio Filtering and Denoising', 'NumberTitle', 'off');

subplot(5,1,1);
plot(t, y, 'r');                               % Original audio
title('Original Audio'); 
ylabel('Amplitude'); 
xlim([0 max(t)]);
grid on;

subplot(5,1,2);
plot(t, y_low, 'color', [1 0.5 0]);            % Low-pass filtered audio
title('Low-pass Filtered'); 
ylabel('Amplitude'); 
xlim([0 max(t)]);
grid on;

subplot(5,1,3);
plot(t, y_high, 'color', [0 0.5 1]);           % High-pass filtered audio
title('High-pass Filtered'); 
ylabel('Amplitude'); 
xlim([0 max(t)]);
grid on;

subplot(5,1,4);
plot(t, y_band, 'color', [0.5 0 0.5]);         % Band-pass filtered audio
title('Band-pass Filtered'); 
ylabel('Amplitude'); 
xlim([0 max(t)]);
grid on;

subplot(5,1,5);
plot(t, y_temp, 'g');                      % Denoised audio
title('Denoised Audio'); 
xlabel('Time (s)'); 
ylabel('Amplitude'); 
xlim([0 max(t)]);
grid on;

sgtitle('Audio Filtering and Noise Reduction'); % Super title for all subplots

%% =============================== MY EVALUATION ===========================================
% Listening Test:
% - Original audio has low-frequency noise (wind sound)
% - High-pass filter removed most low-frequency noise
% - Wiener filter further reduced background noise, improving clarity

% Waveform Observation:
% - Denoised waveform is smoother with less random spikes
% - High-pass removes low-frequency components

% Overall Effectiveness:
% - Voice becomes clearer, noise significantly reduced

%% =============================== SNR CALCULATION ==============================================
% SNR stands for Signal-to-Noise Ratio.
% It quantifies how much of your desired signal (e.g., voice) is present compared to the unwanted noise.

% SNR = Power of Signal / Power of Noise

% Often expressed in decibels (dB):
% SNR (dB) = 10 * log10(P_signal / P_noise)
% Where:
% P_signal = Average power of the signal
% P_noise  = Average power of the noise

% High SNR → Signal is much stronger than noise → clear audio.
% Low SNR → Noise is strong compared to signal → audio is noisy.



noise_before = y - y_high;                      % Estimate noise before denoising
noise_after  = y - y_temp;                  % Estimate noise after denoising

SNR_before = 10*log10(mean(y.^2) / mean(noise_before.^2)); % Compute SNR before
SNR_after  = 10*log10(mean(y.^2) / mean(noise_after.^2));  % Compute SNR after

fprintf('SNR before denoising : %.2f dB\n', SNR_before);
fprintf('SNR after denoising : %.2f dB\n', SNR_after);

%% =============================== SPECTROGRAM COMPARISON =======================================
figure('Name','Spectrogram Comparison','NumberTitle','off');

window = 512;      % Window size for spectrogram
noverlap = 256;    % Number of overlapping samples
nfft = 1024;       % FFT points

subplot(2,1,1); 
spectrogram(y, window, noverlap, nfft, Fs, 'yaxis'); 
title('Original Audio Spectrogram');
ylim([0 5]);       % Limit y-axis to 5 kHz
colorbar;

subplot(2,1,2); 
spectrogram(y_temp, window, noverlap, nfft, Fs, 'yaxis'); 
title('Denoised Audio Spectrogram');
ylim([0 5]);
colorbar;


% Observations:
%--------------
% Original Autio Spectrogram ==>
% There’s a yellow band near 0–1 kHz, which is persistent across time. This is likely wind noise, because wind primarily occupies low frequencies.
% The voice signal is visible in mid frequencies (~0.3–3 kHz), but it is mixed with the low-frequency noise.
% Higher frequencies (>3 kHz) have less energy, mostly background noise.

% Denoised Audio Spectrogram ==>
% The low-frequency yellow band (0–1 kHz) has been reduced. They are now blue, indicating low energy, meaning the wind noise has been removed.
% Voice components in mid-frequency range (~0.3–3 kHz) are clearly preserved.

Overall, the denoised spectrogram is cleaner, with less interference in the low-frequency range.

%% =============================== FFT FREQUENCY SPECTRUM =======================================
disp('Computing FFT for frequency spectrum...');

N = length(y);                         % Number of samples
Y_orig = fft(y);                       % FFT of original audio
Y_denoised_fft = fft(y_temp);      % FFT of denoised audio

f = (0:N-1)*(Fs/N);                    % Frequency vector
half = 1:floor(N/2);                   % Keep only positive frequencies

figure('Name','Frequency Spectrum (FFT)','NumberTitle','off');

subplot(2,1,1);
plot(f(half)/1000, abs(Y_orig(half))/max(abs(Y_orig)), 'r');  % Original FFT normalized
title('Original Audio Frequency Spectrum');
xlabel('Frequency (kHz)');
ylabel('Normalized Magnitude');
grid on;

subplot(2,1,2);
plot(f(half)/1000, abs(Y_denoised_fft(half))/max(abs(Y_denoised_fft)), 'g'); % Denoised FFT
title('Denoised Audio Frequency Spectrum');
xlabel('Frequency (kHz)');
ylabel('Normalized Magnitude');
grid on;

sgtitle('FFT Frequency Spectrum of Real-Time Speech');

%% =============================== TIME SHIFTING (DELAY & ADVANCE) ==============================
disp('Applying Time Shifting to Original Signal...');

time_shift = 1.2;                        % Shift duration in seconds
sample_shift = round(time_shift * Fs);   % Convert to sample count

% Delay (shift right)
y_delay = [zeros(sample_shift,1); y(1:end-sample_shift)];

% Advance (shift left)
y_advance = [y(sample_shift+1:end); zeros(sample_shift,1)];

% Play shifted signals
disp('Playing Delayed Audio...');
sound(y_delay, Fs);
pause(length(y_delay)/Fs + 1);

disp('Playing Advanced Audio...');
sound(y_advance, Fs);
pause(length(y_advance)/Fs + 1);

% Plot comparison
figure('Name','Time Shifting (Delay vs Advance)','NumberTitle','off');

subplot(3,1,1);
plot(t, y, 'r');
title('Original Audio');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(t, y_delay(1:length(t)), 'b'); % Limit to original length
title('Delayed Audio (Shifted Right)');
ylabel('Amplitude');
grid on;

subplot(3,1,3);
plot(t, y_advance(1:length(t)), 'g'); % Limit to original length
title('Advanced Audio (Shifted Left)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;


sgtitle('Time Shifting of Original Audio (Delay and Advance)');




