%% DAAP - Homework A.A. 2024/25
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
disp("================================");
disp("Reading: " + filename);
disp("================================");

% read audio file (and sampling rate?)
[x, fs] = audioread(filename);
x = (x(:,1) + x(:,2))/2;    % from stereo to mono
x = x./max(abs(x));         % normalization
t = [0 : 1/fs : 1/fs*(length(x)-1)].';  % time axis for x(n)

figure();
plot(t,x);
title(filename);

%% Encode LPC coefficients for selected audio file
% work with sampling rate 8 KHz
%fs = 8e3;

% apply a pre-emphasis HPF (freqz(b,1) to plot transfer function)
b = [1, -0.975];    % denominator coefficients
x = filter(b, 1, x);

% LPC parameters
M = 256;
lpc_orders = [4, 10];   % LPC orders (unvoiced and voiced)

% hamming window (length=256, hop=256)
win_len = M;
hop_size = M;
win = hamming(win_len);
n_frames = floor((length(x) - win_len)/hop_size) + 1;

% voiced vs unvoiced frame
[is_voiced, zcr] = voicedframedetection(x, win, hop_size);
p = lpc_orders(is_voiced + 1).';

%% for each frame
for n = 1 : n_frames
    % frame selection and windowing
    frame = x((n-1) * hop_size + 1 : (n-1)*hop_size + win_len);
    frame = frame .* win;

    % compute auto-correlation (size p+1)
    r = xcorr(frame, 'biased');         % auto-correlation function
    r = r((length(r)-1)/2 + 1 : end);   % taking only positive index values
    r = r(1 : p(n));                    % taking the first p values

    % % Wiener-Hopf matrix
    % R = zeros(length(r));
    % for k = 1 : length(r)
    %     for j = 1 : length(r)
    %         R(j,k) = r(abs(j-k));
    %     end
    % end
    % 
    % % LPC coefficients
    % a = R\r;    % a = [R]^-1 * r (size p, same as r)

    % LPC coefficients - alternative
    [a, e, k] = levinson(r, p(n));  % auto-regressive coefficients from
                                    % auto-correlation r and order p(n)
    
    % shaping filter (defined in [0, pi])
    [W, w] = freqz(1, a);

    % prediction error - filter
    e_0 = filter(1, k, frame);
    
end