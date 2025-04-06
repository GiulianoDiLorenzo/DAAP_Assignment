function excite = generateexcitationsignal(voicedIdx, gains, pitch, winLen)
%%generateexcitationsignal(voicedIdx, gains, pitch, winLen)
% This function generates the excitation signal for LPC-10 coding that
% depends on voiced (train of pulses) or unvoiced (random noise) signals
% Arguments:
%   - voicedIdx
%   - gains
%   - pitch
%   - winLen
%
% Outputs
%   - excite
%
% DAAP HW1 2025
% Mirco Pezzoli

    n_frames = length(voicedIdx);
    excite = zeros(n_frames, winLen);
    for n = 1:n_frames
        if voicedIdx(n)
            % For voiced frames: generate a pulse train.
            % The pulse period is given by the pitch period (rounded to an integer number of samples).
            % T = round(pitch(n));
            T = pitch(n);
            pulse = zeros(winLen, 1);
            pulse(1:T:winLen) = 1;  % pulses every T samples
            excite(n, :) = gains(n) * pulse.';
        else
            % For unvoiced frames: generate white noise.
            excite(n, :) = gains(n) * randn(1, winLen);
        end
    end
end
