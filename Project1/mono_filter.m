function Hd = mono_filter
%MONO_FILTER Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.4 and Signal Processing Toolbox 8.0.
% Generated on: 24-Oct-2018 15:28:29

% Chebyshev Type II Lowpass filter designed using FDESIGN.LOWPASS.

% All frequency values are in kHz.
Fs = 256;  % Sampling Frequency

Fpass = 15;          % Passband Frequency
Fstop = 18;          % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 100;         % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
Hd = design(h, 'cheby2', 'MatchExactly', match);

% [EOF]
