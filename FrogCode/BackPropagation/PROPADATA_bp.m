% PROPADATA is a script for propagating retrieved pulses
% through real material.  It plots before and after amplitudes
% and phases in both domains.  A retrieved Speck.Dat file is
% required for the ELOAD command to function.
clc;  clear all
warning off;
PhaseBlank = 0.005;
% 41623VacG.Speck.dat
% Speck42023VacG.Speck.dat
[L2,EL2,~] = ELoad_bp('FilTubeGrenouille_42223fs2_m55000fs3_m3d9e6fs4_Vaccum_85mJ_2019-09-23-175547.speck.dat');
[L,EL,~] = ELoad_bp('Speck42023VacG.Speck.dat');
L0 = L(length(L)./2 +1);
[W,ES] = ToConstFreqSpc_bp(L,EL);
ES = quickscale_bp(ES);
ET = quickscale_bp(ifftc(ES));
dT = 2*pi./(W(end) - W(1));
T = dT.*[-length(W)/2:length(W)/2-1]';
Wavelen = ltow_bp(W);
% % To simulate a typical GRENOUILLE, use:
% direction = -1;
% ESNew = quickscale(betaGREN(ltow(W),ES,direction));
%
% Change the 'direction' to 1 for forward propagation.
% Change the 'direction' to -1 for backward propagation.

% % To simulate just one chunk of material, use:
len = 15;
mat = 'BK7';
ESNew = quickscale_bp(propagate_bp(Wavelen,ES,len,mat));

len = 5.0;
mat = 'FUSED SILICA';
ESNew = quickscale_bp(propagate_bp(Wavelen,ESNew,len,mat));

ETNew = quickscale(ifftc(ESNew));

TemporalFWHMbefore = sprintf('FWHM = %6.1f fs',FWHM_bp(IandP_bp(ET),T));
TemporalFWHMafter = sprintf('FWHM = %6.1f fs', FWHM_bp(IandP_bp(ETNew),T));
SpectralFWHMbefore = sprintf('FWHM = %6.1f nm',FWHM_bp(IandP_bp(ES),ltow_bp(W)));
SpectralFWHMafter = sprintf('FWHM = %6.1f nm', FWHM_bp(IandP_bp(ESNew),ltow_bp(W)));


figure();clf;
% Before (back)propagation
subplot 231;
plotfrog_bp(ET,ET,T,L0);
subplot 232;
plotcmplx_bp(T,ET,[],[],PhaseBlank);
title(TemporalFWHMbefore);
subplot 233
plotcmplx_bp(ltow_bp(W),ES,[],[],PhaseBlank);
% title(SpectralFWHMbefore);

% After (back)propagation
subplot 234;
plotfrog_bp(ETNew,ETNew,T,L0);
subplot 235;
plotcmplx_bp(T,ETNew,[],[],PhaseBlank);
title(TemporalFWHMafter);
subplot 236
plotcmplx(ltow(W),ESNew,[],[],PhaseBlank);
% title(SpectralFWHMafter);

% Colormap HSV2 is for screens, HSV3 is for printing
colormap(hsv3_bp);
warning on;
