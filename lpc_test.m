%% DAAP course 2025 - Homework 1: LPC-based Speak and Spell
% Implement the LPC10 speech encoding used in 1978 Speak and Spell toy by
% Texas Instruments. The TMC0280 chip used to synthesize speech using LPC.
% 
% DAAP HW1 2025
% Mirco Pezzoli

% Di Lorenzo Giuliano
% Longhi Filippo

clc; clear; close all;

addpath('functions')

set(groot,  'DefaultFigureWindowStyle', 'docked', ...       % all figures as tabs in single window
            'DefaultTextInterpreter', 'latex', ...          % interpreter Latex - text and annotations
            'DefaultAxesTickLabelInterpreter', 'latex', ... % interpreter Latex - tick labels
            'DefaultLegendInterpreter', 'latex', ...        % interpreter Latex - legends
            'DefaultLineLineWidth', 1.5, ...                % functions
            'DefaultAxesFontSize', 12, ...                  % axis and title
            'DefaultTextFontSize', 14, ...                  % sgtitle
            'DefaultAxesFontName', 'Times New Roman', ...   % axis and title
            'DefaultTextFontName', 'Times New Roman', ...   % sgtitle
            'DefaultAxesLineWidth', 1, ...                  % axis
            'DefaultConstantLineLineWidth', 1.2, ...        % xline and yline
            'DefaultAxesTitleFontSizeMultiplier', 1.2, ...  % title
            'DefaultFigureColor', 'w', ...                  % background color
            'DefaultAxesBox', 'on', ...                     % plot box
            'DefaultLegendLocation', 'best' ...             % legend position
            );


%% Select an audio file
% Audio file path
filename = input("Choose file name (with right extension): ", "s");

disp("================================");
disp("Reading: " + filename);
disp("================================");

[s, sr] = audioread("input/" + filename);  % read audio file
s = mean(s, 2);                 % from stereo to mono
fs = 8e3;                       % work with sampling rate 8 KHz
s = resample(s, fs, sr);        % resampling original audio
s = s./max(abs(s));             % normalization

%% LPC encoding
encoder;

%% LPC decoding
decoder;

% soundsc(sRec, fs);

% Plot audio file
t = (0 : length(s)-1) / fs;     % time axis for audio file
figure()
plot(t,s)
title(filename)
xlabel("$t$ [s]")
ylabel("$s$")
xlim([min(t) max(t)])
ylim([-1 1])
grid on
