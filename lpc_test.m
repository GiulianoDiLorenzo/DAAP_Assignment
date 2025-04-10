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

[s, sr] = audioread("input/" + filename);  % read audio file
s = mean(s, 2);                 % from stereo to mono
fs = 8e3;                       % work with sampling rate 8 KHz
s = resample(s, fs, sr);        % resampling original audio
s = s./max(abs(s));             % normalization

disp("Reading complete");
disp("================================");

plot_bool = false;      % plots images (use a.mp3)

%% LPC encoding
encoder;

%% LPC decoding
decoder;

%% Comparison and results
% soundsc(s, fs);       % original
% soundsc(s_rec, fs);   % synthesized

% ------------------ PLOTS SECTION ------------------
if plot_bool
    % Plot audio file
    figure(999)
    sgtitle("File " + filename + " - Waveforms and spectra")

    t = (0 : length(s)-1) / fs;
    subplot(2, 2, 1)
    plot(t,s)
    title("Original file - time domain")
    xlabel("$t$ [s]")
    ylabel("$s$")
    xlim([min(t) max(t)])
    ylim([-1 1])
    grid on
        
    N = length(s);
    S = fft(s);
    f = (0:N-1) * (fs / N);
    subplot(2, 2, 2)
    plot(f(1:floor(N/2)), db(abs(S(1:floor(N/2)))))
    title("Original file - frequency domain")
    xlabel("$f$ [Hz]")
    ylabel("$|S|$ [dB]")
    grid on
    
    t = (0:length(s_rec)-1) / fs;
    subplot(2, 2, 3)
    plot(t, s_rec);
    title("Synthesized file - time domain")
    xlabel("$t$ [s]")
    ylabel("$\hat{s}$")
    xlim([min(t) max(t)])
    ylim([-1 1])
    grid on
    
    N = length(s_rec);
    S = fft(s_rec);
    f = (0:N-1) * (fs / N);
    subplot(2, 2, 4)
    plot(f(1:floor(N/2)), db(abs(S(1:floor(N/2)))));
    title("Synthesized file - frequency domain")
    xlabel("$f$ [Hz]")
    ylabel("$|\hat{S}|$ [dB]")
    grid on
end
% ------------------ PLOTS SECTION ------------------