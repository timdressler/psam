% Parameters
fs = 44100;                % Sampling frequency (Hz)
duration = 5;              % Duration of the sine wave (seconds)
freq = 400;                % Frequency of the sine wave (Hz) - A4 note
t = 0:1/fs:duration-1/fs;  % Time vector

% Generate the sine wave for both channels
left_channel = sin(2 * pi * freq * t);   % Left channel sine wave
right_channel = sin(2 * pi * freq * t);  % Right channel sine wave (same frequency)

% Create 1 second of silence
silence = zeros(1, fs);  

% Combine silence, sine wave, and silence again
left_channel = [silence, left_channel, silence];
right_channel = [silence, right_channel, silence];

% Combine the two channels to form a stereo signal
stereo_signal = [left_channel; right_channel]';

% Save the stereo signal as a .wav file
audiowrite([num2str(freq) 'stereo_sine_wave.wav'], stereo_signal, fs);


