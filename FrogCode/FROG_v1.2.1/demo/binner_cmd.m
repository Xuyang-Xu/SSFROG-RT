function Filename = binner_cmd(in_file,SaveDir,sz,width,Method)
    % Command line version of binner GUI. This take in customize parameters to
    % bin the trace. 
    % 'in_file' is the .frg file of the trace
    % 'SaveDir' is the directory to save to (0 for current directory)
    % 'width' : is the bin width in percentage
    % 'sz' : is the size of the grid
    % 'method' : is the method. 1 for fit delay, 2 for fit wavelength

    % By Jeff Wong (GaTech) - 2011-08-09, 2022
    % Edited Fasano - 9/20/2019
    
    % Load the frog trace with the parameters
    [Isig, tau, lam, dtau, lam0, dlam, numDelay, numLam, filename] =...
                    frogload_all(in_file,'delay');

    % Calculate Autocorreatlion and Frequency Convolution as the marginals.
    [ac, freq] = marginals(Isig); 

    % Calculate the Center wavelength and the delay
    lam_c = find_center(lam, freq, 1);
    tau_c = find_center(tau, ac, 3);

    f = ltof(lam);
    tau = tau - tau_c;
    f0 = ltof(lam_c);

    % Binner_CalcGrid
    switch Method
        case 1
            % % case 'Fit Delay'
            dtau = 2 * max(abs([min(tau), max(tau)])) / sz / (width/100);
            tau2 = (-sz/2:sz/2-1) * dtau;
            df = 1 / sz / dtau;
            f2 = (-sz/2:sz/2-1) * df + f0;
            f2 = f2';
            
        case 2
            % % case 'Fit Wavelengths'
            df = (max(f) - min(f)) / sz / (width/100);
            f2 = (-sz/2:sz/2-1) * df + f0;
            f2 = f2';
            dtau = 1 / sz / df;
            tau2 = (-sz/2:sz/2-1) * dtau;
    end

    assignin('base','tau_bincmd',tau);
    assignin('base','lam_bincmd',lam);
    assignin('base','f2_bincmd',f2);
    assignin('base','tau2_bincmd',tau2);

    % Binner_GridData
    % This is the process for binning, basically it is a 2D interpolation
    [Isig, ~, ~] = Binner_GridData(Isig, tau, lam, tau2, f2);
    % Save binned FROG Trace
    if(SaveDir == 0)
        [SaveDir, fname, ~] = fileparts(filename);
    else
        [~, fname, ~] = fileparts(filename);
    end
    Filename = fullfile(SaveDir, sprintf('BIN_%s_b%i_w%i_m%i.bin.frg',fname, sz, width, Method));
    frogsave(quickscale(Isig), dtau, f0, df, [], Filename, Filename);
    fprintf(1,'Finished Binning\n');

end

% Find Center
function y = find_center(X, Y, meth)
    switch meth
        case 1
            y = X(floor(end/2+1));
        case 2
            y = first_moment(Y(:), X(:));
        case 3
            y = X(maxindex(Y(:)));
        otherwise
            error('Unknown centering method.');
    end
end
