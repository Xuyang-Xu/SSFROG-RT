    % Written By Nicholas Fasano
    % 
    % Direct beam - Only one saved trace at 150mJ
    clear all; close all; clc
    
 
    % Assumed that wavelength is the vertical axes (height of image)
    MaxIter = [300];                 % Number of iterations for retrieval
    PlotIter = 1000;                 % iterations to print results
    TemporalCalibration = 2;     % fs/pixel -- corrected value as of 08/22/2023 is 0.455fs/pixel
    SpectralCalibration = 0.4614;     % nm/pixel
    BackGround = 0;                 % order(10) - 8 bit // order(4000) - 16 bit
    CentralWavelength = 400; % 400.16;         % nm
    SZ = 256;                        % size for binner file
    WIDTH = 70;                      % width for binner file

    %========================================= 
    FileDIR =  'Z:\Data\2026_02_15_PlasmaMirrors_PhaseMeasurement\kHz_FROG_ExpTable\'; % Location of traces
    SearchStr = 'duration_kHz_2ndorder_37800fs^2_3rdorder_-54760.0fs^3_20260215_182358.Tiff';
    TifFiles = dir([FileDIR,SearchStr]);
    [~,idx] = sort([TifFiles.datenum],'ascend');
    TifFiles = TifFiles(idx);
    TifFiles = TifFiles(1);
    nshot = length(TifFiles);
    
    % Allocate memory (See Copy_of_GrenouilleAlgorithm_V2.m for definitions)
    Asig = zeros(256,256,nshot); %  
    Esig = zeros(256,256,nshot); %
    tau =  zeros(256,nshot); % time [fs] (x-axis)
    Et =  zeros(256,nshot);  % pulse amplitude in time (a.u.)
    lam1 =  zeros(256,nshot); % frequency [nm] (x-axis)
    lam2 =  zeros(256,nshot);
    w1 =  zeros(256,nshot);
    w2 =  zeros(256,nshot);
    Ew =  zeros(256,nshot); % pulse frequency in time (a.u.)
    
    FWHMt = zeros(nshot,1); % width in time [fs]
    FWHMw = zeros(nshot,1); % width in frequency [nm]
    Energy = zeros(nshot,1);
    ShotNum = zeros(nshot,1);
    Z = zeros(MaxIter,nshot);
    G = zeros(MaxIter,nshot);
%     
    
for j = 1:nshot
        [Asig(:,:,j), Esig(:,:,j), tau(:,j), Et(:,j), lam1(:,j), lam2(:,j), w1(:,j),...
            w2(:,j), Ew(:,j), FWHMt(j), FWHMw(j), Z(:,j), G(:,j)] =...
                        Copy_of_GrenouilleAlgorithm_V2(MaxIter, PlotIter, TemporalCalibration,...
                        SpectralCalibration,BackGround,CentralWavelength,SZ,WIDTH,...
                        TifFiles(j));
end    



    % ==============================================
    % ==============================================
    % ========== Post data calculations ============
    % ==============================================
    % ==============================================
    
    It = zeros(SZ,nshot); Iw = It;
    phiet = It;     phiew = It;
    etmax = zeros(nshot,1);   jetmax = etmax;
    ewmax = etmax;          jewmax = etmax;
    for j = 1:nshot
        It(:,j) = magsq(Et(:,j)); Iw(:,j) = magsq(Ew(:,j));
        [etmax(j), jetmax(j)] = max(It(:,j));
        phiet(:,j) = -unwrap(angle(Et(:,j)));  

        [ewmax(j), jewmax(j)] = max(Iw(:,j));
        phiew(:,j) = -unwrap(angle(Ew(:,j)));
    end
    

     % ==============================================
    % ==============================================
    % === Plot Original and reconstructed traces ===
    % ==============================================
    % ==============================================     
    figure('Position',[24 123 1319 504])
    for j = 1:nshot
        clf
        subplot(1,3,1)
        imagesc(tau(:,j), lam2(:,j), Asig(:,:,j))
        colormap(frogcolormap)
        caxis([0,1])
        %ylim([341, 459]); %xlim([-250, 250])
        set(gca,'YDir','normal')
        title('Original')
        xlabel('Delay [fs]'); ylabel('Wavelength [nm]')

        subplot(1,3,2)
        imagesc(tau(:,j), lam2(:,j), abs(Esig(:,:,j)))
        colormap(frogcolormap)
        caxis([0,1])
        %ylim([341, 459]); xlim([-250, 250])
        set(gca,'YDir','normal')
        title('Reconstructed')
        xlabel('Delay [fs]');% ylabel('Frequency')

        subplot(1,3,3)
        imagesc(tau(:,j), lam2(:,j), abs(Esig(:,:,j)) - Asig(:,:,j))
        colormap(frogcolormap)
        colorbar
        %ylim([341, 459]); xlim([-250, 250])
        set(gca,'YDir','normal')
        title('Difference')
        xlabel('Delay [fs]');% ylabel('Frequency')
        
        pause(0.1)

    end

%%
    % ==============================================
    % ==============================================
    % ========== Plot pulse and spectrum ============
    % ==============================================
    % ==============================================
    CMAP = [0.0000 0.4470 0.7410;
            0.8500 0.3250 0.0980];
    hfigOne = figure('Position',[6 91 651 576]);
    for j = 1:nshot
            clf;
            
            subplot(2,1,1)
            yyaxis left
            h1 = plot(tau(:,j),(It(:,j)/etmax(j)),'-','Color',CMAP(1,:));
            hold on
            ylabel('Intensity'); xlabel('Delay [fs]')  
            xlim([-150,150]); %ylim([.9e-5,1])

            yyaxis right
            h2 = plot(tau(:,j),phiet(:,j) - phiet(jetmax(j),j),'--','Color',CMAP(2,:));
            hold on           
            ylabel('Phase [rad]');
            xlim([-150,150]); %ylim([-1,1])
            text(-140,.75,['FWHM = ',num2str(FWHMt(j),'%.2f'),'fs'],'FontSize',14,'FontWeight','Bold')
            
            subplot(2,1,2)
            yyaxis left
            h1 = plot(lam1(:,j),Iw(:,j)/ewmax(j),'-','Color',CMAP(1,:));
            hold on       
            xlabel('Frequency [nm]'); ylabel('Intensity'); %ylim([.9e-5,1])
            xlim([725, 875]);
            
            yyaxis right
            h2 = plot(lam1(:,j),phiew(:,j) - phiew(jewmax(j),j),'--','Color',CMAP(2,:));
            hold on     
            xlim([725, 875]); %ylim([-1,1])     
            ylabel('Phase [rad]'); xlabel('Frequency [nm]')
            text(730,.75,['FWHM = ',num2str(FWHMw(j),'%.2f'),'nm'],'FontSize',14,'FontWeight','Bold')
    
            pause(0.1)
    
    end

% 
% %% 
% % FileDIR =  'Z:\Data\2025_05_20_Oscillator_FROG\FROG_Traces\'; % Location of traces
% %    SearchStr = '1_Oscillator_0.2um_1024_20250520.tiff';
% wavelength_crop3 = lam1(89:147);
% spectral_phase3 = phiew(:,j) - phiew(jewmax(j),j);
% spectral_phase_crop3 = spectral_phase3(89:147);
% quadratic3 = fittedmodel3(wavelength_crop3);
% % x_val = linspace(wavelength_crop2(end), wavelength_crop2(1), length(wavelength_crop2));
% % quadratic2 = zeros(1, length(wavelength_crop2));
% % for i=1:length(wavelength_crop2)
% %     quadratic2(i) = 0.01861*x_val(i)^2-5.29*x_val(i);
% % end
% difference3 = spectral_phase_crop3-flip(quadratic3);
% figure; plot(wavelength_crop3, difference3)
% xlabel('Wavelength [nm]'); ylabel('Spectral phase [rad]')
% %title('Spectral phase - (c*x^2+d*x) fit')
% box on
% set(gca, 'FontSize', 20)
% 
% figure; plot(wavelength_crop3,spectral_phase_crop3); hold on; plot(wavelength_crop3, flip(quadratic3)); legend('Exp',  'Fit')
% 
% % val = 
% % 
% %      General model:
% %      val(x) = a*x^4+b*x^3+c*x^2+d*x
% %      Coefficients (with 95% confidence bounds):
% %        a =   8.452e-09  (7.659e-09, 9.245e-09)
% %        b =  -2.175e-05  (-2.369e-05, -1.982e-05)
% %        c =     0.01861  (0.01703, 0.02018)
% %        d =       -5.29  (-5.716, -4.864)
% 

