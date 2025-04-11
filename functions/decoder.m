% LPC-decoder
% 
% DAAP course 2025
% Mirco Pezzoli
%
% The following script synthesize an audio file
% starting from the LPC-10 encoding that includes
% ['fs', 'win_len', 'hop_size', 'n_frames', 'lpc_coeffs', 'gains', 'pitch_periods', 'is_voiced']

%% Generate excitation signal
% Generate train of impulses (voiced) or white noise (unvoiced)
u = generateexcitationsignal(is_voiced, gains, pitch_periods, win_len);


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
   
    % Shaping filtering H
    frame_synth = filter(1, A, u(n, :).');
    
    % Place synthesized frame into right output position
    s_synth((n-1)*hop_size + 1 : (n-1)*hop_size + win_len) = frame_synth;


    % ------------------ PLOTS SECTION ------------------
    if plot_bool
        if ismember(n, plot_idx)
            % Plots time frame and spectra - FIGURE N+1
            figure(n+1)
            
            % Frame in time
            t = (0 : length(frame_synth)-1) / fs;
            subplot(2, 2, 3)
            plot(t, frame_synth)
            title("Reconstructed frame - time domain")
            xlabel("$t$ [ms]")
            ylabel("$\hat{s}_n$")
            grid on
            xlim([min(t) max(t)])
    
            % Frame in frequency
            S = fft(frame_synth);
            f = (0:win_len-1)*(fs/win_len);
            subplot(2, 2, 4)
            plot(f(1:win_len/2), db(abs(S(1:win_len/2))))
            title("Reconstructed frame - frequency domain")
            xlabel("$f$ [Hz]")
            ylabel("$|\hat{S}_n|$ [dB]")
            grid on
            xlim([0 fs/2])
    
            % Plot excitation signals - FIGURE N+2
            figure(n+2) 
            subplot(2,1,2);
            plot(t, u(n, :));
            title("Excitation signal")
            xlabel("$t$ [ms]")
            ylabel("$e'$")
            grid on
            xlim([min(t) max(t)])
        end
    end
    % ------------------ PLOTS SECTION ------------------
end

% A de-emphasis filter is applied to reverse the pre-emphasis done at the encoder.
% b_deemp = abs([1, -0.975]);
s_rec = filter(abs(b), 1, s_synth);
s_rec = s_rec ./ max(s_rec);

disp("Decoding complete");
disp("================================");