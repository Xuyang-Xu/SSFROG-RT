    % Written By Nicholas Fasano
    % 
    % Direct beam - Only one saved trace at 150mJ
    clear all; close all; clc
    
 
    % Assumed that wavelength is the vertical axes (height of image)
    MaxIter = [100];                  % Number of iterations for retrieval
    PlotIter = 1000;                 % iterations to print results
    TemporalCalibration = 0.455;     % fs/pixel -- corrected value as of 08/22/2023 is 0.455fs/pixel
    SpectralCalibration = 0.055;     % nm/pixel
    BackGround = 0;                  % order(10) - 8 bit // order(4000) - 16 bit
    CentralWavelength = 390.3;       % nm
    SZ = 256;                        % size for binner file
    WIDTH = 70;                      % width for binner file

    %========================================= 
    FileDIR =  'Z:\Data\2024_06_17_FROG\GRENOUILLE_10Hz\'; % Location of traces
    SearchStr = '10Hz_GRENOUILLE_7*.tiff';
    TifFiles = dir([FileDIR,SearchStr]);
    [~,idx] = sort([TifFiles.datenum],'ascend');
    TifFiles = TifFiles(idx);
    TifFiles = TifFiles(4); %Was 10
    nshot = length(TifFiles);
    
    % Allocate memory (See GrenouilleAlgorithm_V2.m for definitions)
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
    
    
for j = 1:nshot
        [Asig(:,:,j), Esig(:,:,j), tau(:,j), Et(:,j), lam1(:,j), lam2(:,j), w1(:,j),...
            w2(:,j), Ew(:,j), FWHMt(j), FWHMw(j), Z(:,j), G(:,j)] =...
                        GrenouilleAlgorithm_V2(MaxIter, PlotIter, TemporalCalibration,...
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
            xlim([-150,150]); ylim([.9e-5,1])

            yyaxis right
            h2 = plot(tau(:,j),phiet(:,j) - phiet(jetmax(j),j),'--','Color',CMAP(2,:));
            hold on           
            ylabel('Phase [rad]');
            xlim([-150,150]); ylim([-1,1])
            text(-140,.75,['FWHM = ',num2str(FWHMt(j),'%.2f'),'fs'],'FontSize',14,'FontWeight','Bold')
            
            subplot(2,1,2)
            yyaxis left
            h1 = plot(lam1(:,j),Iw(:,j)/ewmax(j),'-','Color',CMAP(1,:));
            hold on       
            xlabel('Frequency [nm]'); ylabel('Intensity'); ylim([.9e-5,1])
            xlim([725, 875]);
            
            yyaxis right
            h2 = plot(lam1(:,j),phiew(:,j) - phiew(jewmax(j),j),'--','Color',CMAP(2,:));
            hold on     
            xlim([725, 875]); ylim([-1,1])     
            ylabel('Phase [rad]'); xlabel('Frequency [nm]')
            text(730,.75,['FWHM = ',num2str(FWHMw(j),'%.2f'),'nm'],'FontSize',14,'FontWeight','Bold')
    
            pause(0.1)
    
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
        ylim([380, 480]); xlim([-250, 250])
        set(gca,'YDir','normal')
        title('Original')
        xlabel('Delay [fs]'); ylabel('Wavelength [nm]')

        subplot(1,3,2)
        imagesc(tau(:,j), lam2(:,j), abs(Esig(:,:,j)))
        colormap(frogcolormap)
        caxis([0,1])
        ylim([380, 480]); xlim([-250, 250])
        set(gca,'YDir','normal')
        title('Reconstructed')
        xlabel('Delay [fs]');% ylabel('Frequency')

        subplot(1,3,3)
        imagesc(tau(:,j), lam2(:,j), abs(Esig(:,:,j)) - Asig(:,:,j))
        colormap(frogcolormap)
        colorbar
        ylim([380, 480]); xlim([-250, 250])
        set(gca,'YDir','normal')
        title('Difference')
        xlabel('Delay [fs]');% ylabel('Frequency')
        
        pause(0.1)

    end








