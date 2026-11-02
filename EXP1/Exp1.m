%% EXPERIMENT 1 - ECG NOISE REMOVAL AND QRS DETECTION
clear;
close all;
clc;

%% INPUT PARAMETERS
fs = 360;  % Sampling frequency in Hz
filename = 'ecg_sample.mat';

%% CHECK AND LOAD THE ECG FILE
if ~isfile(filename)
    error(['File not found: ', filename, ...
        '. Place the file in the MATLAB Current Folder.']);
end

S = load(filename);

if isfield(S, 'ecg')
    ecg = S.ecg;
elseif isfield(S, 'ECG')
    ecg = S.ECG;
else
    error('ECG variable not found. Check the MAT file contents.');
end

% Convert ECG to a column vector
ecg = double(ecg(:));

% Check ECG data
if isempty(ecg) || any(~isfinite(ecg))
    error('ECG data is empty or contains invalid values.');
end

if length(ecg) < 100
    error('ECG signal is too short for reliable analysis.');
end

t = (0:length(ecg)-1)' / fs;

fprintf('Samples = %d\n', length(ecg));
fprintf('Duration = %.2f seconds\n', t(end));

%% FIGURE 1: RAW ECG
figure('Name', 'Raw ECG');
plot(t, ecg, 'k');
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Raw ECG Signal');

%% PREPROCESSING
ecg0 = ecg - mean(ecg);

%% BAND-PASS FILTER
lowCut = 5;
highCut = 15;
filterOrder = 4;

[b, a] = butter(filterOrder, ...
    [lowCut highCut] / (fs/2), 'bandpass');

ecgFiltered = filtfilt(b, a, ecg0);

%% FIGURE 2: FILTERED ECG
figure('Name', 'Filtered ECG');
plot(t, ecgFiltered, 'b');
grid on;
xlabel('Time (s)');
ylabel('Filtered amplitude');
title('ECG after Band-Pass Filtering');

%% QRS ENHANCEMENT
dECG = [0; diff(ecgFiltered)];
squared = dECG .^ 2;

integrationWindow = round(0.150 * fs);
integrated = movmean(squared, integrationWindow);

%% FIGURE 3: QRS ENHANCEMENT
figure('Name', 'QRS Enhancement');
plot(t, integrated);
grid on;
xlabel('Time (s)');
ylabel('Integrated energy');
title('Moving-Window Integrated ECG');

%% PEAK DETECTION
threshold = 0.35 * max(integrated);
minimumDistance = round(0.25 * fs);

if max(integrated) <= 0
    qrsLocs = [];
else
    [~, qrsLocs] = findpeaks(integrated, ...
        'MinPeakHeight', threshold, ...
        'MinPeakDistance', minimumDistance);
end

%% FIGURE 4: QRS DETECTION
figure('Name', 'Detected QRS');
plot(t, ecgFiltered, 'b');
hold on;

if ~isempty(qrsLocs)
    plot(t(qrsLocs), ecgFiltered(qrsLocs), 'ro', ...
        'MarkerFaceColor', 'r');
end

grid on;
xlabel('Time (s)');
ylabel('ECG amplitude');
title('Filtered ECG with Detected QRS Candidates');
legend('Filtered ECG', 'QRS candidates');

%% HEART RATE CALCULATION
if length(qrsLocs) >= 2
    RR = diff(qrsLocs) / fs;
    meanRR = mean(RR);
    HR = 60 / meanRR;
else
    RR = [];
    meanRR = NaN;
    HR = NaN;

    warning(['Fewer than two QRS candidates detected. ', ...
        'Heart rate cannot be estimated reliably.']);
end

%% DISPLAY RESULTS
fprintf('\n===== OBSERVATIONS =====\n');
fprintf('Sampling frequency = %.2f Hz\n', fs);
fprintf('Filter range = %.2f to %.2f Hz\n', lowCut, highCut);
fprintf('Detected QRS candidates = %d\n', length(qrsLocs));

if isnan(meanRR)
    fprintf('Mean RR interval = Not available\n');
    fprintf('Estimated heart rate = Not available\n');
else
    fprintf('Mean RR interval = %.4f s\n', meanRR);
    fprintf('Estimated heart rate = %.2f bpm\n', HR);
end
