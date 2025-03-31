function pitch = pitchdetectionamdf(e)
%pitchdetectionamdf
% This function find the pitch of a frame, using the Average Magnitude
% Difference Function AMDF.
% Arguments
%   - e
%
% Output
%   - pitch
%
% Based on AMDF see https://www.researchgate.net/publication/
% 228854783_Pitch_detection_algorithm_autocorrelation_method_and_AMDF
%
% DAAP HW1 2025
% Mirco Pezzoli

winLen = length(e);
len =length(e);
lags = winLen;
amd = zeros(1,lags);
e = [e; zeros(lags,1)];
for k=1:lags
    for j=1:len
        % see equation 6 of reference
        amd(k) = 
    end
    amd(k) = 
end
pitch = find(amd == min(amd(25:80)));
end


