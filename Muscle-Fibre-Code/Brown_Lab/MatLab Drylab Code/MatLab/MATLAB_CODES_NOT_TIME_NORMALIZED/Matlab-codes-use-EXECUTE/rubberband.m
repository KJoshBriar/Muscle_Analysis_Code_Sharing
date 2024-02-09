function [N] = rubberband(nFr, M)
% N - rubberbanded output data.
% nFr - the number of frames in the rubberbanded output time-series.
% M - the original time series data (can have multiple columns).


[R,C] = size(M);
N = zeros(R,C);
nFr = R; %Bypass rubberband after the fact with Global outlength;
rat = (R-1)/(nFr-1);
for i = 1 : nFr
    s = 1 + rat*(i-1);
    if i == nFr 
        j = floor(s-0.5);
    else    
        j = round(s-0.5);
    end    
    f = s-j;
    N(i,:) = (1-f)*M(j,:) + f*M(j+1,:);
end  
end
