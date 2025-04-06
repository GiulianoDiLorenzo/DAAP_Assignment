% LPC-decoder
% 
% DAAP course 2025
% Mirco Pezzoli
%
% The following script synthesize an audio file
% starting from the LPC-10 encoding that includes
% ['filename_short', 'fs', 'M', 'n_frames', 'a', 'lpc_coeffs', 'gains', 'pitch_periods', 'is_voiced']

%% Load the encoded data
% What about passing all necessary info inside .mat?
% Like fs, M, n_frames, 

% Load the encoded parameters from the encoder:
% load("lpc10_" + filename_short + ".mat");

% Also, define the sampling rate (fs) as used in encoding (e.g., 8000 Hz)
% fs = 8000;
% 
% M = 256;            % Frame length (same as in encoder)
% n_frames = size(lpc_coeffs, 1);  % Number of frames

% lpc_orders = [4, 10];   % LPC orders (unvoiced and voiced)

% %[
% is_voiced(isnan(is_voiced)) = 0;
% gains(isnan(gains)) = 0;
% pitch_periods(isnan(pitch_periods)) = 0;
% lpc_coeffs(isnan(lpc_coeffs)) = 0;
% %]

%% Generate excitation signal
% Generate train of impulses (voiced) or white noise (unvoiced)
excite = generateexcitationsignal(is_voiced, gains, pitch_periods, M);

% Plot excitation signals
figure; 
subplot(2,1,1);
plot(excite(42, :)); % frame 80
title('Excitation Signal for Frame 80');
subplot(2,1,2);
plot(excite(80, :)); % frame 42
title('Excitation Signal for Frame 42');


% Pre-allocate the reconstructed signal.
synthS = zeros(n_frames * M, 1);

%% Decode audio signal
disp("================================");
disp("Decoding: " + filename_short + ".mp3");

for n = 1:n_frames
    p = lpc_orders(is_voiced(n) + 1);
    % % Determine LPC order for this frame:
    % if is_voiced(n)
    %     p = 10;  % Voiced frames were modeled with order 10
    % else
    %     p = 4;   % Unvoiced frames were modeled with order 4
    % end    
    
    % Construct the LPC synthesis filter A(z) = 1 - a1*z^(-1) - a2*z^(-2) - ... - ap*z^(-p)
    A = [1, -lpc_coeffs(n, 1:p)];
    
    % Retrieve the excitation signal for this frame:
    u = excite(n, :).';  % ensure it is a column vector of length M
    
    % Synthesize the frame by filtering the excitation through the LPC filter.
    frame_synth = filter(1, A, u);
    
    % Place the synthesized frame into the output signal.
    synthS((n-1)*M + 1 : (n-1)*M + M) = frame_synth;
    
    % Plot only frame 40 and frame 84:
    if ismember(n, [42, 80])
        figure;
        t = (0:1/fs:(M-1)/fs);  % time axis for one frame
        plot(t, frame_synth);
        title(sprintf('Synthesized Frame %d', n));
        xlabel('Time (s)');
        ylabel('Amplitude');
        grid on;
    end
end

% A de-emphasis filter is applied to reverse the pre-emphasis done at the encoder.
b_deemp = abs([1, -0.975]);
sRec = filter(b_deemp, 1, synthS);

% Create a time axis for the reconstructed signal:
t = (0:length(sRec)-1) / fs;  % fs is the sampling frequency

% Plot the time-domain signal:
figure;
plot(t, synthS);
xlabel('Time (s)');
ylabel('Amplitude');
title('Reconstructed Speech Signal');

% Optionally, plot the frequency spectrum:
N = length(sRec);
S = fft(sRec);
f = (0:N-1) * (fs / N);
figure;
plot(f(1:floor(N/2)), abs(S(1:floor(N/2))));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum of Reconstructed Signal');



disp("Decoding complete");
disp("================================");