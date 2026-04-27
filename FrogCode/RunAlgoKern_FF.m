function [Et, Esig, Z, G] = RunAlgoKern_FF(Asig,Et,Esig)
    % ==== Nicholas Fasano
    % ==== Last Edited: 09/19/2019
    % ==== Descrption: Main algorithm for FROG calculation using second
    % ==== harmonic generation approach. This is a modified version from Trebino's
    % ==== website (http://frog.gatech.edu/).

    % If the FrogObj.Esig hasn't been initialized, do so.
    if(Esig == 0)
        Esig = CalcEsig(Et, Et);
% %         Esig = fft_FROG(Esig);
        Esig = fftc(Esig, [], 1);
    end

    % The magnitude replacement step.
    Esig = MagRepl_Edited(Esig, Asig);

    % Take Esig to the time and delay domain.
    Esig = ifftc(Esig,[],1);

    % Get the new guess and Z error for E(t).
    dZ = -dZdE_shg(Esig,Et);
    [Et, Z] = MinZerr_shg(Esig, Et, dZ);
    Et = center(Et,'max');

    % Generate a new Esig.
    Esig = CalcEsig(Et, Et);

    % Take Esig back to the frequency and delay domain.
    Esig = fftc(Esig, [], 1);

    % Calculate the G error and correction factor.
    [G,a] = MinGerr(Esig, Asig);

    % Apply the correction factor to E(t).
    Et = Et .* a^(1/6);
end
