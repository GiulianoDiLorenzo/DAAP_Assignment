function [voicedIdx, zc] = voicedframedetection(s,win,hopsize)
%% voicedframedetection(s,win,hopsize)
% Estimate if a frame is voiced or unvoiced based on zerocrossing rate.
% Arguments
%   - s
%   - win
%   - hopsize
%
% Outputs
%   - voicedIdx
%   - zc
%
% DAAP HW1 2025
% Mirco Pezzoli

winLen = length(win);
sLen = length(s);
nFrame = floor( (sLen - winLen)/hopsize ) + 1;
voicedIdx = zeros(nFrame, 1);
zc = voicedIdx;

for ii =1:nFrame
    fIdx = (ii-1)*hopsize+1 : (ii-1)*hopsize+winLen;
    sn = s(fIdx).*win;
    for j =2:winLen
        temp(j-1) = abs(sign(sn(j))-sign(sn(j-1)));
    end
    zc(ii)=sum(temp);
end
for ii = 1:nFrame
    fIdx = (ii-1)*hopsize+1 : (ii-1)*hopsize+winLen;
    sn = s(fIdx).*win;
    if (zc(ii)>mean(zc))                      % Detecting Zero crossing
        voicedIdx(ii)=0;
    elseif  (sum(sn(:).^2)<0.001)
        voicedIdx(ii) = 0;
    else
        voicedIdx(ii)=1;
    end
end