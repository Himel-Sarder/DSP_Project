clc;
clear all;
close all;

disp('Downloading your voice file...');
url = 'https://raw.githubusercontent.com/Himel-Sarder/DSP_Project/main/lab.wav';
filename = 'lab.wav';
websave(filename, url);
disp('Download complete!');

% Read audio
[y, Fs] = audioread(filename);


% Play original audio
disp("Playing Original Voice..")
sound(y, Fs);

% Time vector
t = (0:length(y)-1)/Fs;

% Plot waveform
figure;
plot(t, y, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('Waveform of Original Voice');

%% Low-pass Filter
disp('Applying Low-pass Filter...');
Fc_low = 1000;
[b_low, a_low] = butter(6, Fc_low/(Fs/2), 'low');
y_low = filter(b_low, a_low, y);
disp('Playing Low-pass Filtered Voice...');
sound(y_low, Fs);
pause(length(y)/Fs + 1);

%% High-pass Filter
disp('Applying High-pass Filter...');
Fc_high = 1000;
[b_high, a_high] = butter(6, Fc_high/(Fs/2), 'high');
y_high = filter(b_high, a_high, y);
disp('Playing High-pass Filtered Voice...');
sound(y_high, Fs);
pause(length(y)/Fs + 1);

%% Band-pass Filter
disp('Applying Band-pass Filter...');
Fc1 = 500; Fc2 = 2000;
[b_band, a_band] = butter(6, [Fc1 Fc2]/(Fs/2), 'bandpass');
y_band = filter(b_band, a_band, y);
disp('Playing Band-pass Filtered Voice...');
sound(y_band, Fs);
pause(length(y)/Fs + 1);

%% Plot all signals
figure;
subplot(4,1,1); 
plot(t, y, 'r'); 
title('Original Audio'); 
xlabel('Time (s)'); 
ylabel('Amplitude');


subplot(4,1,2); 
plot(t, y_low); 
title('Low-pass Filtered'); 
xlabel('Time (s)'); 
ylabel('Amplitude');


subplot(4,1,3); 
plot(t, y_high); 
title('High-pass Filtered'); 
xlabel('Time (s)'); 
ylabel('Amplitude');


subplot(4,1,4); 
plot(t, y_band); 
title('Band-pass Filtered'); 
xlabel('Time (s)'); 
ylabel('Amplitude');

