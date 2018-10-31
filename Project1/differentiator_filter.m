function Hd = differentiator_filter
%DIFFERENTIATOR_FILTER Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.4 and DSP System Toolbox 9.6.
% Generated on: 16-Oct-2018 08:46:08

% Equiripple Differentiator filter designed using the FIRPM function.

% All frequency values are normalized to 1.

N = 100;                       % Order
F = [0 0.95 0.98 1];  % Frequency Vector
A = [0 1 0 0];                 % Amplitude Vector
W = [1 1];                     % Weight Vector

% Calculate the coefficients using the FIRPM function.
b  = firls(N, F, A, W, 'differentiator');
Hd = dfilt.dffir(b);
[gd, ] = grpdelay(Hd);
disp(gd(1));
%fvtool(Hd)
% [EOF]
