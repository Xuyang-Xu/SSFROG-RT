function [T, ET, ETNew, Wavelen, ES, ESNew, FWHMTo, FWHMTf, FWHMWo, FWHMWf] =...
    BackPropagation(SpeckFileName,LenBK7,LenFusedSi)

    % Nicholas Fasano
    % Modified PropData function for adjusting retrieved pulse to account for
    % dispersion through glass.
    % if either LenBK7 or LenFusedSi < 0, then iterate lengths until
    % shortest pulse is found --> Bisection methos is used here.
    % Here an effective Bk7 Length is found and LenFusedSi = 0.
    
    % EL is spectral field as a function of wavelength
    % ES = spectral field as a function of frequency
    % ET and ETNEW are the fields in the time domain
    % Note: L and Wavelen are equivalent -- wavelength in nm of main beam
    
    [Wavelen,EL,~] = ELoad_bp(SpeckFileName);
    [W,ES] = ToConstFreqSpc_bp(Wavelen,EL);
    ET = ifftc(ES);
    dT = 2*pi./(W(end) - W(1));
    T = dT.*[-length(W)/2:length(W)/2-1]';
    

% % % % %     if(LenBK7 < 0 || LenFusedSi < 0)
% % % % %         % iterate LenBK7 until optimal pulse length is found
% % % % %         len_guess = linspace(0,100,100);
% % % % %         FWHM_guess = zeros(100,1);   
% % % % %         material = 'BK7';
% % % % %         for j = 1:length(len_guess)
% % % % %             ESNew = propagate_bp(Wavelen,ES,len_guess(j),material);
% % % % %             ETNew = ifftc(ESNew);
% % % % %             FWHM_guess(j) = FWHM_bp(IandP_bp(ETNew),T);
% % % % %         end
% % % % %         figure, plot(len_guess,FWHM_guess)
% % % % %         hold on
% % % % %         xlabel('BK7 Length [mm]'); ylabel('FWHM [fs]')
% % % % %         [~,jfwhm] = min(FWHM_guess);
% % % % %         ESNew = propagate_bp(Wavelen,ES,len_guess(jfwhm),material);
% % % % %         ETNew = ifftc(ESNew);
% % % % %         fprintf('LenBK7 for Ideal pulse = %.2fmm\n',len_guess(jfwhm));
% % % % %     else
        % used supplied lengths
        
        % Adjust for BK7 glass 
        len = LenBK7; % length [mm]
        material = 'BK7';
        ESNew = propagate_bp(Wavelen,ES,len,material);
        
% %         % Adjust for Fused Silica glass 
% %         len = LenFusedSi; % length [mm]
% %         material = 'FUSED SILICA';
% %         ESNew = propagate_bp(Wavelen,ESNew,len,material);
        ETNew = ifftc(ESNew);
% % % % %     end
    
    FWHMTo = FWHM_bp(IandP_bp(ET),T);
    FWHMTf = FWHM_bp(IandP_bp(ETNew),T);
    FWHMWo = FWHM_bp(IandP_bp(ES),Wavelen);
    FWHMWf = FWHM_bp(IandP_bp(ESNew),Wavelen);
end
