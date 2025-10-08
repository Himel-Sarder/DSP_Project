%% Complete Signal Analysis Script for NoisyVoice.wav

clc; clear; close all;

%% Step 1: Load the audio file
filename = 'NoisyVoice.wav';
[y, Fs] = audioread(filename);  % y = audio samples, Fs = sampling rate
if size(y,2) > 1
    y = mean(y,2); % convert stereo to mono
end
N = length(y);       % number of samples
t = (0:N-1)/Fs;      % time vector
duration = N/Fs;     % duration in seconds

%% Step 2: Audio file info
info = audioinfo(filename);
disp('--- Audio File Info ---');
disp(info);

%% Step 3: Basic statistics
fprintf('\n--- Signal Statistics ---\n');
fprintf('Number of samples: %d\n', N);
fprintf('Duration: %.2f seconds\n', duration);
fprintf('Mean amplitude: %.4f\n', mean(y));
fprintf('Variance: %.4f\n', var(y));
fprintf('Standard deviation: %.4f\n', std(y));
fprintf('Max amplitude: %.4f\n', max(y));
fprintf('Min amplitude: %.4f\n', min(y));

%% Step 4: Time-domain plot
figure;
plot(t, y);
xlabel('Time (s)');
ylabel('Amplitude');
title('Time-domain Signal');
grid on;

%% Step 5: Listen to the signal
disp('Playing original audio...');
sound(y, Fs);

%% Step 6: Frequency-domain analysis (FFT)
Y = fft(y);
f = (0:N-1)*(Fs/N);   % frequency vector
magnitude = abs(Y)/N;

figure;
plot(f, magnitude);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Spectrum');
xlim([0 Fs/2]);       % Nyquist limit
grid on;

%% Step 7: Spectrogram (Time-Frequency)
figure;
spectrogram(y, 256, 250, 256, Fs, 'yaxis');
title('Spectrogram');
