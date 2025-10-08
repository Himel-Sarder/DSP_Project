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

% Convert to mono if stereo
if size(y,2) == 2
    y = mean(y, 2);  % average left and right channels
end

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

%% ---------------------- STUDENT EVALUATION ----------------------
% Listening Test:
% The original audio had noticeable low-frequency noise, like wind sound.
% After applying the high-pass filter, most of the low-frequency noise was removed,
% but some high-frequency noise and small artifacts remained.
% Using the Wiener filter on top of the high-pass filtered signal further reduced
% the background noise, making the voice clearer and more intelligible.

% Waveform Observation:
% Comparing the waveforms, the denoised audio shows smoother amplitude variations
% and less random noise spikes compared to the original.
% The high-pass filtered waveform removed low-frequency components,
% while the Wiener filter smoothed out smaller fluctuations effectively.

% Overall Effectiveness:
% The combination of high-pass filtering and Wiener denoising worked well for this audio sample.
% The voice is clearer, and the background noise is significantly reduced.

%% ---------------------- SNR CALCULATION ----------------------
noise_before = y - y_high;      
noise_after  = y - y_denoised; 

SNR_before = 10*log10(mean(y.^2) / mean(noise_before.^2));
SNR_after  = 10*log10(mean(y.^2) / mean(noise_after.^2));

fprintf('SNR before denoising : %.2f dB\n', SNR_before);
fprintf('SNR after denoising : %.2f dB\n', SNR_after);

%% ---------------------- SPECTROGRAM COMPARISON ----------------------
figure('Name','Spectrogram Comparison','NumberTitle','off');

window = 512;      
noverlap = 256;    
nfft = 1024;       

subplot(2,1,1); 
spectrogram(y, window, noverlap, nfft, Fs, 'yaxis'); 
title('Original Audio Spectrogram');
ylim([0 5]);  
colorbar;

subplot(2,1,2); 
spectrogram(y_denoised, window, noverlap, nfft, Fs, 'yaxis'); 
title('Denoised Audio Spectrogram');
ylim([0 5]);
colorbar;

%% ---------------------- FFT FREQUENCY SPECTRUM ----------------------
disp('Computing FFT for frequency spectrum...');

% FFT for Original and Denoised signals
N = length(y);
Y_orig = fft(y);
Y_denoised_fft = fft(y_denoised);

% Frequency axis
f = (0:N-1)*(Fs/N);
half = 1:floor(N/2);  

figure('Name','Frequency Spectrum (FFT)','NumberTitle','off');

subplot(2,1,1);
plot(f(half)/1000, abs(Y_orig(half))/max(abs(Y_orig)), 'r');
title('Original Audio Frequency Spectrum');
xlabel('Frequency (kHz)');
ylabel('Normalized Magnitude');
grid on;

subplot(2,1,2);
plot(f(half)/1000, abs(Y_denoised_fft(half))/max(abs(Y_denoised_fft)), 'g');
title('Denoised Audio Frequency Spectrum');
xlabel('Frequency (kHz)');
ylabel('Normalized Magnitude');
grid on;

sgtitle('FFT Frequency Spectrum of Real-Time Speech');

%% ---------------------- MANUAL DFT (Will Take Huge Time) ----------------------
% disp('Computing Manual DFT for frequency spectrum...');

% Function for manual DFT
% function X = manual_dft(x)
%     N = length(x); 
%     X = zeros(1, N); 
    
%    for k = 0:N-1
%        sum_val = 0;
%        for n = 0:N-1
            % W_N = exp(-j * 2 * pi / N)
%             exponent = -1j * 2 * pi * n * k / N;
%             sum_val = sum_val + x(n+1) * exp(exponent); 
%         end
%         X(k+1) = sum_val;
%     end
% end

% N = length(y); 
% Y_orig_dft = manual_dft(y');
% Y_denoised_dft = manual_dft(y_denoised'); 

% Frequency axis
% f = (0:N-1)*(Fs/N); 
% half = 1:floor(N/2); 

