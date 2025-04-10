% LPC-encode
% 
% DAAP course 2025
% Mirco Pezzoli
% 
% The following script computes the LPC coefficients
% for each frame of an audio file with sr=8000

%% Define parameters
% Apply a pre-emphasis HPF (freqz(b,1) to plot transfer function)
b = [1, -0.975];        % denominator coefficients
s = filter(b, 1, s);

% Hanning windowing (window length 256 samples, 50% overlap)
win_len = 256;      
hop_size = win_len/2;
win = hanning(win_len);
n_frames = floor((length(s) - win_len)/hop_size) + 1;

% LPC model parameters
lpc_coeffs = zeros(n_frames, 10);
pitch_periods = zeros(n_frames, 1);
gains = zeros(n_frames, 1);
frame_errors = zeros(n_frames, win_len);

% Voiced vs unvoiced frames
[is_voiced, ~] = voicedframedetection(s, win, hop_size);

% Custom plots
plot_idx = [42, 82];

%% Estimate LPC coefficients
disp("================================");
disp("Encoding: " + filename);

for n = 1 : n_frames
    % Frame selection and windowing
    frame = s((n-1)*hop_size + 1 : (n-1)*hop_size + win_len) .* win;

    % LPC order
    if is_voiced(n)
        p = 10;
    else
        p = 4;
    end

    % Compute autocorrelation
    r = xcorr(frame, 'biased');  
    r = r(win_len:end);  % only positive lags
    r = r(1 : p+1);   % r(0) and first p values

    % A(z) coefficients using Levinson-Durbin recursion
    [A, ~, ~] = levinson(r, p);

    % Take care of NaNs (because of r=0 for silent frames)
    if any(isnan(A))
        % LPC coefficients not updated
        A = zeros(1, p);

    else
        % LPC coefficients, excluding 1
        lpc_coeffs(n, 1:p) = -A(2:end);

        % Prediction error - whitening filtering
        frame_errors(n, :) = filter(A, 1, frame);       

        % voiced vs unvoiced frame
        if is_voiced(n)
            % Low-pass at 800 Hz for pitched frames
            frame_errors(n, :) = lowpass(frame_errors(n, :), 800, fs);

            % Pitch and gain computation after LPF
            pitch_periods(n) = pitchdetectionamdf(frame_errors(n, :).');

            lim = floor(win_len./pitch_periods(n)).*pitch_periods(n);
            power(n) = (1/lim).*(frame_errors(n, 1:lim) * frame_errors(n, 1:lim).');
            gains(n) = sqrt(power(n)*pitch_periods(n));
            
        else
            % MSE of frame as excitation gain
            gains(n) = sqrt(mean(frame_errors(n, :) .^ 2));
            % gains(n) = r(1) - lpc_coeffs(n, 1:p)*r(2:end);
            % pred_gain(n) = r(1) / gains(n);
        end
    end
    
    % ------------------ PLOTS SECTION ------------------
    if plot_bool
        if ismember(n, plot_idx)
            % Plots for spectra, error in time and frequency - FIGURE "N"
            [H, w] = freqz(1, A);   % shaping filter and w axis
            S = fft(frame);         % frame spectrum
            f = (0:win_len-1)*(fs/win_len);     % frequency axis (up to fs)
    
            figure(n)
            sgtitle("File " + filename + " - Frame n." + n + " (p = " + p +")")
    
            % Spectra comparison
            subplot(3,1,1)
            plot(f(1:win_len/2), db(abs(S(1:win_len/2))))
            hold on
            plot(w./(2*pi) * fs, db(abs(H)))
            hold off
            title("Original spectrum $S$ and LPC approximation $H$")
            xlabel("$f$ [Hz]")
            ylabel("$|S_n|$, $|H|$ [dB]")
            grid on
            xlim([0 fs/2])
            legend("Original spectrum", "Shaping filter approximation")
    
            % Error in time
            t = (0 : 1/fs : 1/fs*(win_len-1))*10e3;   % [ms] time axis
            subplot(3,1,2)
            plot(t, frame_errors(n, :))
            title("Prediction error - time domain")
            xlabel("$t$ [ms]")
            ylabel("$e$")
            grid on
            xlim([min(t) max(t)])
    
            % Error in frequency
            E = fft(frame_errors(n, :));    % error spectrum
            subplot(3,1,3)
            plot(f(1:win_len/2), db(abs(E(1:win_len/2))))
            title("Prediction error - frequency domain")
            xlabel("$f$ [Hz]")
            ylabel("$|E|$ [dB]")
            grid on
            xlim([0 fs/2])
        
            % Plots time frame and spectra - FIGURE N+1
            figure(n+1)
            sgtitle("File " + filename + " - Frame n." + n)
    
            % Frame in time
            subplot(2, 2, 1)
            plot(t, frame)
            title("Original frame - time domain")
            xlabel("$t$ [ms]")
            ylabel("$s_n$")
            grid on
            xlim([min(t) max(t)])
    
            % Frame in frequency
            subplot(2, 2, 2)
            plot(f(1:win_len/2), db(abs(S(1:win_len/2))))
            title("Original frame - frequency domain")
            xlabel("$f$ [Hz]")
            ylabel("$|S_n|$ [dB]")
            grid on
            xlim([0 fs/2])
    
            % Plots prediction error and excitation (decoder) - FIGURE N+2
            figure(n+2)
            sgtitle("File " + filename + " - Frame n." + n)
    
            % Error in time
            t = (0 : 1/fs : 1/fs*(win_len-1))*10e3;   % [ms] time axis
            subplot(2,1,1)
            plot(t, frame_errors(n, :))
            title("Prediction error - time domain")
            xlabel("$t$ [ms]")
            ylabel("$e$")
            grid on
            xlim([min(t) max(t)])
        end
    end
    % ------------------ PLOTS SECTION ------------------
end

%% Save encoded data
save("lpc10_encoding.mat", 'fs', 'win_len', 'hop_size', 'n_frames', 'lpc_coeffs', 'gains', 'pitch_periods', 'is_voiced');
disp("Encoding complete");
disp("================================");