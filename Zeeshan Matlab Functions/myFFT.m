% this function takes a vector of the variable (in time) 'x' and time step 
% 'dt' to transform in Foruier domain using MATLAB fft function. It returns
% the frequency signal and the variable in frequency domain.

function [f, X] = myFFT(x,dt)

fs = 1/dt;              
XX = fft(x);
ind = 1:length(XX)/2+1;
X = abs(XX(ind)/fs);
f = [0:1:fs/2];
end

