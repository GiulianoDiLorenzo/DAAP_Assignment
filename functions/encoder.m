%% LPC-encode
%
% DAAP course 2025
% Mirco Pezzoli
% clear
close all
clc

%% Define parameters
% apply a pre-emphasis HPF (freqz(b,1) to plot transfer function)
b = [1, -0.975];    % denominator coefficients
x = filter(b, 1, x);

% hamming window (length=256, hop=256)
M = 256;                % segment length (windows)
win_len = M;
hop_size = M;
win = hamming(win_len);
n_frames = floor((length(x) - win_len)/hop_size) + 1;

% LPC model parameters
lpc_orders = [4, 10];   % LPC orders (unvoiced and voiced)
lpc_coeffs = zeros(n_frames, lpc_orders(2));
pitch_periods = zeros(n_frames, 1);
gains = zeros(n_frames, 1);
frame_errors = zeros(n_frames, win_len);

% voiced vs unvoiced frames
[is_voiced, zcr] = voicedframedetection(x, win, hop_size);
p = lpc_orders(is_voiced + 1).';

%% Estimate LPC coefficients
disp("================================");
disp("Encoding: " + filename);

for n = 1 : n_frames
    % frame selection and windowing
    frame = x((n-1)*hop_size + 1 : (n-1)*hop_size + win_len) .* win;
    
    % compute autocorrelation
    r = xcorr(frame, 'biased');  
    r = r(win_len:end);  % only positive lags
    r = r(1 : p(n)+1);   % r(0) and first p values

    % LPC coefficients using Levinson-Durbin recursion
    [a, ~, ~] = levinson(r, p(n));
    lpc_coeffs(n, 1:p(n)) = a(2:end);   % LPC coefficients, excluding 1
    
    % ----------- PEZZOLI -----------
    % Alternative solution use lpc
    % a(ii, 1:orderLPC+1) = lpc(sn.*win, orderLPC);

    % Avoid NaNs
    % ----------- PEZZOLI -----------

    % prediction error - whitening filtering
    frame_errors(n, :) = filter(a, 1, frame);

    % MSE of frame - excitation gain
    gains(n) = sqrt(mean(frame_errors(n, :) .^ 2));
    
    % pitch period computation for voiced frames
    if (is_voiced(n))
        pitch_periods(n) = pitchdetectionamdf(frame_errors(n, :));
    end
    
    % ----------- PEZZOLI -----------
    % Plot the Magnitude of the signal spectrum, and on the same graph, the
    % the LPC spectrum estimation (remember the definition of |E(omega)|)
    if doPlot
        
        % Compute the shaping filter H using freqz function
         
        %frequency axis
        % Compute the DFT of the original signal
        %FFT of the signal
        

        figure(3), clf
        subplot(3,1,1)
        %Plot of the FFT of the original signal
        hold on;
        title();
        xlabel();
        %Spectral matching of the filter
        hold off

        subplot(3,1,2), hold on
       
       % plot prediction error (time domain)
        xlabel('Time [ms]');
        title()
        % Plot the prediction error magnitude spectrum
        subplot(3,1,3)
        %Plot prediction error (frequency domain)
        title ('');
        
    end
    % Pitch estimation for voiced signals
    if voicedIdx(n) == 1
        
        % Optionally low pass the error
        
    end

    % Plot the predition error e in time
    if voicedIdx(n) == 1
        % Compute the gain for a pitched sound
        lim = floor(winLen./pitch(n)).*pitch(n);
        power(n) = (1/lim).*(e(1:lim)'*e(1:lim));
        gains(n) = sqrt(power(n)*pitch(n));

    else
        % Compute the gain for unvoiced sounds
        
    end
    % ----------- PEZZOLI -----------
end

%% Save encoded data
save('lpc10_encoded.mat', 'lpc_coeffs', 'gains', 'pitch_periods', 'is_voiced');
disp("Encoding complete");
disp("================================");