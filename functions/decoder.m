%% LPC decoder

% clear 
% close all
% clc

%% Load the encoded data

% Load the encoded parameters from the encoder:
load('lpc10_encoded.mat', 'lpc_coeffs', 'gains', 'pitch_periods', 'is_voiced');

% Also, define the sampling rate (fs) as used in encoding (e.g., 8000 Hz)
fs = 8000;

M = 256;            % Frame length (same as in encoder)
n_frames = size(lpc_coeffs, 1);  % Number of frames


%[
is_voiced(isnan(is_voiced)) = 0;
gains(isnan(gains)) = 0;
pitch_periods(isnan(pitch_periods)) = 0;
lpc_coeffs(isnan(lpc_coeffs)) = 0;
%]

% Generate excitation signals for each frame.
% This function should create a pulse train for voiced frames (using pitch_periods)
% and white noise for unvoiced frames, both scaled by the corresponding gain.
excite = generateexcitationsignal(is_voiced, gains, pitch_periods, M);

figure; 
subplot(2,1,1);
plot(excite(42, :)); % frame 80
title('Excitation Signal for Frame 80');
subplot(2,1,2);
plot(excite(80, :)); % frame 42
title('Excitation Signal for Frame 42');


% Pre-allocate the reconstructed signal.
synthS = zeros(n_frames * M, 1);

disp("================================");
disp("Decoding: ");

for n = 1:n_frames
    % Determine LPC order for this frame:
    if is_voiced(n)
        p = 10;  % Voiced frames were modeled with order 10
    else
        p = 4;   % Unvoiced frames were modeled with order 4
    end
    
    % Construct the LPC synthesis filter A(z) = 1 - a1*z^(-1) - a2*z^(-2) - ... - ap*z^(-p)
    A = [1, -lpc_coeffs(n, 1:p)];
    
    % Retrieve the excitation signal for this frame:
    u = excite(n, :).';  % ensure it is a column vector of length M
    
    % Synthesize the frame by filtering the excitation through the LPC filter.
    frame_synth = filter(1, A, u);
    
    % Place the synthesized frame into the output signal.
    start_idx = (n-1)*M + 1;
    synthS(start_idx:start_idx+M-1) = frame_synth;
    
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
    %print also in encoder to check!!!
end




% A de-emphasis filter is applied to reverse the pre-emphasis done at the encoder.
b_deemp = abs([1, -0.975]);
sRec = filter(b_deemp, 1, synthS);

soundsc(sRec, fs)

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