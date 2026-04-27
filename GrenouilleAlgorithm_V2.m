%=======================================================%
%============ Nicholas Fasano ==========================%
%============ Analyze Grenouille Trace =================%
%============ Last Modified: 09/24/2019 ================%
%============ Adapted from Rick Trebino ================%
%============ http://frog.gatech.edu/code.html==========%
%=======================================================%

function [Asig, Esig, tau, Et, lam1, lam2, w1, w2, Ew, FWHMt, FWHMw, Z, G] =...
                GrenouilleAlgorithm_V2(MaxIter, PlotIter, TemporalCalibration,...
                                       SpectralCalibration,BackGround,...
                                       CentralWavelength,SZ,WIDTH,TraceFileInfo)

    % Add Frog code to working matlab path
    FrogDir = '\\lockhart.princeton.edu\elmil\Code\Experimental_Analysis\FROG\Frog_Automated\FrogCode';
    addpath(genpath(FrogDir));
% INPUTS:
%     MaxIter                  % Number of iterations for retrieval
%     PlotIter                 % iterations to print results
%     TemporalCalibration      % fs/pixel
%     SpectralCalibration      % nm/pixel
%     BackGround               % Additional background (an automatic
                               % background subtraction is performed within
                               % GrenouilleAlgorithm_V2.m) 
%     CentralWavelength        % nm
%     SZ                       % size for binner file
%     WIDTH                    % width for binner file
%     TraceFileInfo            % returned structure from dir() of file info   

% OUTPUS:
    % Asig                     % Original Binned Trace
    % Esig                     % Retrieved Binned Trace
    % tau                      % Delay (temporal axis) [fs]
    % Et, Ew                   % Temporal, spectral fields (complex numbers)
    % lam1, lam2               % Wavelength axis [nm] for fundamental,second harmonic
    % w1, w2                   % frequency axis [nm] for fundamental,second harmonic
    % FWHMt, FWHMw             % Temporal [fs], frequency [nm] FWHM 
    % Z, G                     % Two Error measures indicating successful
                               % convergence of the retrieval algorithm
    
    
    for j = 1

        % ==============================================
        % ==============================================
        % ========== Load in original trace, copy ======
        % ========== it locally, and do background =====
        % ========== subtraction and filtering =========
        % ==============================================
        % ==============================================  
        TraceFileName = TraceFileInfo.name;    
        LocalDir = ['BinnedFiles\'];

        if(~exist(LocalDir,'dir')), mkdir(LocalDir); end

        Source = [TraceFileInfo.folder,'\',TraceFileName];
        Destination = [LocalDir,'\',TraceFileName];
        copyfile(Source, Destination);  

        OriginalImage = imread(Destination);
        OriginalImage = OriginalImage(2:end,2:end);
        
        % Find suitable background subtraction
        % Method One: Average of noise from a single row far from the trace
        SortBackground = sort(OriginalImage(30,:));
        MeanBackground = mean(SortBackground(end-500:end));
        
        BinnedImage = OriginalImage - MeanBackground - BackGround; 
        BinnedImage = flipud((double(medfilt2(BinnedImage,[5 5]))));
        BinnedImage(BinnedImage < 0) = 0;

        % ==============================================
        % ==============================================
        % ========== Create .frg file in delay format====
        % ==============================================
        % ==============================================

        % get Wavlength and temporal width from image
        NumPixelsWavelength = length(BinnedImage(:,1));
        NumPixelsTemporal   = length(BinnedImage(1,:));

        % Create .frg file with same name as .tif
        header1 = [NumPixelsTemporal, NumPixelsWavelength,...
            TemporalCalibration, SpectralCalibration, CentralWavelength];
        outfile = [LocalDir,'\',TraceFileName(1:end-4),'.frg'];
        dlmwrite(outfile,header1,'delimiter','\t')
        dlmwrite(outfile,BinnedImage,'delimiter','\t','-append')

        % ==============================================
        % ==============================================
        % ========== Bin .frg file with binner GUI =====
        % ========== or binner_cmd function ============
        % ==============================================
        % ==============================================

        % Method One -- Call binner GUI
        % binner

        % Method Two -- Use binner_cmd
        BinnedFileName =  binner_cmd(outfile,LocalDir,SZ,WIDTH,1);

        % ==============================================
        % ==============================================
        % ========== retrieve pulse with frogger GUI ===
        % ========== or FastFrogger function ===========
        % ==============================================
        % ==============================================

        % Method One -- Call froger GUI
        % frogger

        % Method 2: Use FastFrogger    
        [Asig, Esig, tau, Et, lam1, lam2, w1, w2, Ew, FWHMt, FWHMw, Z, G] = ...
                    FastFrogger(BinnedFileName, MaxIter, PlotIter);


    end
end
