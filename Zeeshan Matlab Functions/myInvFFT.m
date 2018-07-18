% This function performs inverse Fourier transform of the harmonic
% coefficients given in input X. 

% input -->  a) X contains static and Fourier coefficients of cosines only
%            b) nStep is arbitrary time discretization parameter% 
% output --> x is the vector in time domain. 
%            another FFT of the corresponding quantities requires
%            multiplication by 2 because of fft function definition in
%            MATLAB

function x = myInvFFT(X, nStep)

nH = length(X)-1;       % calculate number of coefficients of cosines
XX = zeros(nStep,1);    % define a vector in frequency domain
XX(1) = X(1);           % assign the static value

% place the cosine coefficients at their position
XX(2:length(X)) = X(2:length(X))/2;    
% also place the coefficients at aliasing frequencies
XX(nStep+1-nH:nStep) = conj(X(length(X):-1:2))/2;

x = ifft(nStep*XX);
end