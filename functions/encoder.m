%% LPC-encode
%
% DAAP course 2025
% Mirco Pezzoli
% clear
close all
clc

%% Define parameters
     % Reading input audio signal
% Work at 8kHz sampling rate

% Normalize the signal with respect to peak value


% Apply a pre-emphasis filter as highpass filter
b = [1, -0.975];
s= filter(b,1,s);                   % pre-emphasis filtering

% STFT Parameters
% suggested parameters 256-sample window length and hopsize, hamming window

% LPC-10 parameters

% Pitch estimation parameters
pitch = zeros(nFrame, 1);
% Optionally low pass the prediction error with FIR filter with cut off at
% 800 Hz

doPlot = false; % Activate plots


%% Esimate voiced frame


%% Estimate LPC filters

for ii = 1:nFrame
    % Extract a frame from the signal
    
    % Compute autocorrelation
    
    % Compute the whitening filter parameters using the function levinson
    % (see doc levinson)
    
    % Alternative solution use lpc
    % a(ii, 1:orderLPC+1) = lpc(sn.*win, orderLPC);

    % Avoid NaNs
    
    % Compute the prediction error using filter
    
    

    % Compute MSE of the frame
    

    
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
    if voicedIdx(ii) == 1
        
        % Optionally low pass the error
        
    end

    % Plot the predition error e in time
    if voicedIdx(ii) == 1
        % Compute the gain for a pitched sound
        lim = floor(winLen./pitch(ii)).*pitch(ii);
        power(ii) = (1/lim).*(e(1:lim)'*e(1:lim));
        gains(ii) = sqrt(power(ii)*pitch(ii));

    else
        % Compute the gain for unvoiced sounds
        
    end
end

%% Save encoded data
save()