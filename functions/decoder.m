% LPC-decoder
% 
% DAAP course 2025
% Mirco Pezzoli
%
% The following script synthesize an audio file
% starting from the LPC-10 encoding that includes
% ['fs', 'M', 'n_frames', 'a', 'lpc_coeffs', 'gains', 'pitch_periods', 'is_voiced']

%% Load the encoded data
% % Load the encoded parameters:
% load("lpc10_" + filename_short + ".mat");

%% Generate excitation signal
% Generate train of impulses (voiced) or white noise (unvoiced)
u = generateexcitationsignal(is_voiced, gains, pitch_periods, M);

% Plot excitation signals
figure; 
subplot(2,1,1);
plot(t, u(42, :)); % frame 80
title('Excitation Signal for Frame 42');
xlabel("$t$ [ms]")
subplot(2,1,2);
plot(t, u(80, :)); % frame 42
title('Excitation Signal for Frame 80');
xlabel("$t$ [ms]")

% Pre-allocate the reconstructed signal.
s_synth = zeros(n_frames * M, 1);

%% Decode audio signal
disp("================================");
disp("Decoding: " + filename);

for n = 1:n_frames
    % LPC order
    if is_voiced(n)
        p = 10;
    else
        p = 4;
    end    
    
    % Whitening filter coefficients
    A = [1, -lpc_coeffs(n, 1:p)];
    
    % Excitation for frame n
    % u(n, :).'
    % column vector of length M
    
    % Shaping filtering H
    frame_synth = filter(1, A, u(n, :).');
    
    % Place synthesized frame into right output position
    s_synth((n-1)*M + 1 : (n-1)*M + M) = frame_synth;

    % Plots time frame and spectra
    if ismember(n, plot_idx)
        figure()
        sgtitle("File " + filename + " - Reconstructed frame n." + n)

        % Frame in time
        t = (0 : length(frame_synth)-1) / fs;     % time axis for audio file
        subplot(2, 1, 1)
        plot(t, frame_synth)
        xlabel("$t$ [ms]")
        ylabel("$\hat{s}_n$")
        grid on
        xlim([min(t) max(t)])

        % Frame in frequency
        S = fft(frame_synth);   % frame spectrum
        f = (0:M-1)*(fs/M);     % frequency axis (up to fs)
        subplot(2, 1, 2)
        plot(f(1:M/2), db(abs(S(1:M/2))))
        xlabel("$f$ [Hz]")
        ylabel("$|\hat{S}_n|$ [dB]")
        grid on
        xlim([0 fs/2])
    end
    
    % % Plot only frame 40 and frame 84:
    % if ismember(n, [42, 80])
    %     figure;
    %     t = (0:1/fs:(M-1)/fs);  % time axis for one frame
    %     plot(t, frame_synth);
    %     title(sprintf('Synthesized Frame %d', n));
    %     xlabel('Time (s)');
    %     ylabel('Amplitude');
    %     grid on;
    % end
end

% A de-emphasis filter is applied to reverse the pre-emphasis done at the encoder.
% b_deemp = abs([1, -0.975]);
s_rec = filter(abs(b), 1, s_synth);
s_rec = s_rec ./ max(s_rec);

% Create a time axis for the reconstructed signal:
t = (0:length(s_rec)-1) / fs;  % fs is the sampling frequency

% Plot the time-domain signal:
figure;
plot(t, s_rec);
xlabel('Time (s)');
ylabel('Amplitude');
title('Reconstructed Speech Signal');

% Optionally, plot the frequency spectrum:
N = length(s_rec);
S = fft(s_rec);
f = (0:N-1) * (fs / N);
figure;
plot(f(1:floor(N/2)), abs(S(1:floor(N/2))));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum of Reconstructed Signal');



disp("Decoding complete");
disp("================================");