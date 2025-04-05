% LPC-encode
% 
% DAAP course 2025
% Mirco Pezzoli
% 
% The following script computes the LPC coefficients
% for each frame of an audio file with sr=8000

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

% % customized plots
% doPlot = zeros(n_frames, 1);
% doPlot(80) = 1;
% doPlot(42) = 1;

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
        pitch_periods(n) = pitchdetectionamdf(frame_errors(n, :).');

        % % ----------- PEZZOLI -----------
        % Optionally low pass the error
        % % ----------- PEZZOLI -----------
    end

    % % ----------- PEZZOLI -----------
    % if voicedIdx(n) == 1
    %     % Compute the gain for a pitched sound
    %     lim = floor(winLen./pitch(n)).*pitch(n);
    %     power(n) = (1/lim).*(e(1:lim)'*e(1:lim));
    %     gains(n) = sqrt(power(n)*pitch(n));
    % 
    % else
    %     % Compute the gain for unvoiced sounds
    % 
    % end
    % % ----------- PEZZOLI -----------
    
    % plots for spectra, error in time and frequency
    if doPlot(n)
        
        [H, w] = freqz(1, a);   % shaping filter and w axis
        S = fft(frame);         % frame spectrum
        f = (0:M-1)*(fs/M);     % frequency axis (up to fs)

        figure()
        sgtitle("File " + filename_short + ".mp3 - Frame n." + n + " (is voiced = " + is_voiced(n) +")")

        % spectra comparison
        subplot(3,1,1)
        plot(f(1:M/2), db(abs(S(1:M/2))))
        hold on
        plot(w./(2*pi) * fs, db(abs(H)))
        hold off
        title("Original spectrum $S$ and LPC approximation $H$")
        xlabel("$f$ [Hz]")
        ylabel("$|S|$, $|H|$ [dB]")
        grid on
        xlim([0 fs/2])
        legend("Original spectrum", "Shaping filter approximation")

        % error in time
        t = (0 : 1/fs : 1/fs*(M-1))*10e3;   % [ms] time axis
        subplot(3,1,2)
        plot(t, frame_errors(n, :))
        title("Prediction error - time domain")
        xlabel("$t$ [ms]")
        ylabel("$e$")
        grid on
        xlim([min(t) max(t)])

        % error in frequency
        E = fft(frame_errors(n, :));    % error spectrum
        subplot(3,1,3)
        plot(f(1:M/2), db(abs(E(1:M/2))))
        title("Prediction error - frequency domain")
        xlabel("$f$ [Hz]")
        ylabel("$|E|$ [dB]")
        grid on
        xlim([0 fs/2])
        
    end
end

%% Save encoded data
save("lpc10_" + filename_short + ".mat", 'lpc_coeffs', 'gains', 'pitch_periods', 'is_voiced');
disp("Encoding complete");
disp("================================");