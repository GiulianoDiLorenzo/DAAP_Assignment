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
% audio file path
filename = "input/" + input("Choose file name (no extension): ", "s") + ".mp3";
filename_short = erase(filename, ["input/", ".mp3"]);

disp("================================");
disp("Reading: " + filename);
disp("================================");

[x, sr] = audioread(filename);  % read audio file
x = mean(x, 2);                 % from stereo to mono
fs = 8e3;                       % work with sampling rate 8 KHz
x = resample(x, fs, sr);        % resampling original audio
x = x./max(abs(x));             % normalization

t = [0 : 1/fs : 1/fs*(length(x)-1)].';  % time axis for audio file
figure();
plot(t,x);
title(filename);

%% LPC encoding
encoder;

%% Choose the audio to be decoded
decoder;


