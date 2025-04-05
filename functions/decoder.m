%% LPC decoder

clear 
close all
clc

%% Load the encoded data

% Load the encoded parameters from the encoder:
load('lpc10_encoded.mat', 'lpc_coeffs', 'gains', 'pitch_periods', 'is_voiced');

% Also, define the sampling rate (fs) as used in encoding (e.g., 8000 Hz)
fs = 8000;

M = 256;            % Frame length (same as in encoder)
n_frames = size(lpc_coeffs, 1);  % Number of frames

% Generate excitation signals for each frame.
% This function should create a pulse train for voiced frames (using pitch_periods)
% and white noise for unvoiced frames, both scaled by the corresponding gain.
excite = generateexcitationsignal(is_voiced, gains, pitch_periods, M);



% Pre-allocate the reconstructed signal.
synthS = zeros(n_frames * M, 1);

for n = 1:n_frames
    % Determine LPC order for this frame:
    if is_voiced(n)
        p = 10;  % Voiced frames were modeled with order 10
    else
        p = 4;   % Unvoiced frames were modeled with order 4
    end
    
    % Construct the LPC synthesis filter A(z) = 1 - a1*z^(-1) - a2*z^(-2) - ... - ap*z^(-p)
    % (Our stored coefficients are the a_k's; we prepend 1 and change the sign)
    A = [1, -lpc_coeffs(n, 1:p)];
    
    % Retrieve the excitation signal for this frame:
    u = excite(n, :).';  % ensure it is a column vector of length M
    
    % Synthesize the frame by filtering the excitation through the LPC filter.
    % This implements: s_hat(n) = filter(1, A(z), u)
    frame_synth = filter(1, A, u);
    
    % Place the synthesized frame into the output signal.
    % Here, we assume no overlap (hop_size = frame length).
    start_idx = (n-1)*M + 1;
    synthS(start_idx:start_idx+M-1) = frame_synth;
end




% A de-emphasis filter is applied to reverse the pre-emphasis done at the encoder.
b_deemp = abs([1, -0.975]);
sRec = filter(b_deemp, 1, synthS);

soundsc(sRec, fs)
