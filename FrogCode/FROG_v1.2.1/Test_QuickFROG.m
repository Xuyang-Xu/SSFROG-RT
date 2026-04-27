function Test_QuickFROG(varargin)
% Testing QuickFROG
% This file is written to test the Matlab environment to make sure all the
% necessary files to run the FROG algorithm are present. Please follow the
% readme file in the home directory if errors show up in test run.

close all; clear all; clc

N = 128;
dt = 4;

% [N, dt, domain, plse, err] = parsevarargin(varargin, N, dt, domain, plse, err);

% Pulse Generation 
[Et0, t, Ew0, w] = pulsegenerator(N, @fgaussian, 20, dt, 800,[],0,[]);

Et0 = center(Et0,'max');
Et0 = Et0/Et0(end/2);
Ew0 = fftc(Et0);

% Calculate FROG Trace
Et0 = ifftc(Ew0);
Esig = CalcEsig(Et0,Et0);
Asig = abs(fft_FROG(Esig));
Asig = quickscale(Asig);

% Add noise
Asig = nonegatives(Asig);

Et = pulsegenerator(N); Et = complex(abs(Et));

% Run the Retrieval Algorithm
Et = QuickFROG_tT(Asig, Et, t, w);

FWHM = fwhm(magsq(Et),t);
end
