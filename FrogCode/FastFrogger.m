function [Asig, Esig, tau, Et, lam1, lam2, w1, w2, Ew, FWHMt, FWHMw, Z, G] = ...
                    FastFrogger(BinnedFileName, MaxIter, PlotIter)
    % ==== Nicholas Fasano
    % ==== Last Edited: 09/19/2019
    % ==== Descrption: Main algorithm for FROG calculation using second
    % ==== harmonic generation approach. This is a modified version from Trebino's
    % ==== website (http://frog.gatech.edu/).

    % Constants
    c = 299792458;
    
    % Initializations
    Esig = 0;
    Z = zeros(MaxIter,1);
    G = zeros(MaxIter,1);

    % Load in binned file and initial trace
    [Asig,tau,freq,dtau,f0,df,NumD,NumL,filename] = frogload(BinnedFileName); %NumD = 256, NumL = 256
    
    % We only want the amplitude.
    Asig = sqrt(Asig);

    % HACK!!!!
    % Binnined traces were not saved in expected units.  This hack fixes it.
    w2 = (freq) * 2 * pi;
    lam2 = 2*pi*c./(w2*1e15) * 1e9;
    WaveLen =  ltow(w2 - w2(floor(end/2)+1) / 2); 
    w1 = 2*pi*c./(WaveLen*1e-9);
    lam1 = WaveLen;  

  
    
    % Set a default new first guess.
    Et = frand(tau, tau(floor(end*3/4)));
%     Et = dlmread('Eto.txt');
    % Run one step to get the proper error values.
    [ Et, Esig,Z(1),G(1)] = RunAlgoKern_FF(Asig,Et,Esig);
    freq = w2/(2*pi);
    Ew = fftc(Et);      
% % %     % Plot Initial guess and spectrum
% % %     hh = figure('Position',[309 107 757 559]);
% % %     subplot(2,1,1)
% % %     htime = plotcmplx(tau,Et);
% % %     htime(1).YLabel.Position(1) = htime(1).YLabel.Position(1) - 20;
% % %     htime(2).YLabel.Position(1) = htime(2).YLabel.Position(1) + 20;
% % %     htime(1).Title.Position(2) = htime(1).Title.Position(2) +.02;
% % %     title('Temporal Domain')
% % %     hold off
% % % 
% % %     subplot(2,1,2)
% % %     hfreq = plotcmplx(freq,Ew);
% % %     hfreq(1).YLabel.Position(1) = hfreq(1).YLabel.Position(1) - 10/1000;
% % %     hfreq(2).YLabel.Position(1) = hfreq(2).YLabel.Position(1) + 5/1000;
% % %     hfreq(1).Title.Position(2) = hfreq(1).Title.Position(2) +.02;
% % %     hfreq(1).XLim = [.7,.78];
% % %     hfreq(2).XLim = [.7,.78];
% % %     hfreq(2).YLim = [-10,10];
% % %     figure(hh);
% % %     title('Frequency Domain')  
% % %     hold off
    tic
    if(PlotIter < MaxIter)
        hhh = figure('Position',[6 91 651 576]);
    end
    for j = 2:MaxIter     
        if((mod(j,PlotIter) == 0 || j == 2) && PlotIter < MaxIter)
            % Perform some calcs on time pulse
            It = magsq(Et); Iw = magsq(Ew);
            [etmax, jetmax] = max(It);
            phiet = -unwrap(angle(Et));  

            [ewmax, jewmax] = max(Iw);
            phiew = -unwrap(angle(Ew));          
            Ew = fftc(Et); 
            FWHMt = fwhm(magsq(Et),tau);
            FWHMw = fwhm(magsq(Ew),WaveLen);
              
            % ==============================================
            % ==============================================
            % ========== Plot pulse and spectrum ============
            % ==============================================
            % ==============================================
            
            figure(hhh)
            clf
            subplot(2,1,1)
            yyaxis left
            plot(tau,It/etmax,'-','Color',[0,0,1])
            hold on
            ylabel('Amplitude'); xlabel('Delay [fs]')  
            xlim([-150,150]); ylim([.9e-5,1])

            title(sprintf('Iteration = %d',j))

            yyaxis right
                semilogy(tau,phiet - phiet(jetmax),'--','Color',[0,0,1])
                hold on           
            ylabel('Phase [rad]');
            xlim([-200,200]); ylim([-2,2])

            subplot(2,1,2)
            yyaxis left
                plot(lam1,Iw/ewmax,'-','Color',[0,0,1])
                hold on       
            xlabel('Frequency [nm]'); ylabel('Amplitude'); ylim([.9e-5,1])
            xlim([725, 875]);

            yyaxis right
            semilogy(lam1,phiew - phiew(jewmax),'--','Color',[0,0,1])
            hold on     
            xlim([725, 875]); ylim([-2,2])     
            ylabel('Phase [rad]'); xlabel('Frequency [nm]') 
            

            fprintf('Iteration = %d\t\tFWHM_t = %.2ffs\t\tFWHM_w = %.2fnm\t\tTime_sec = %.2f\n',...
            j,FWHMt,FWHMw,toc) 
        end    

        [Et, Esig, Z(j), G(j)] = RunAlgoKern_FF(Asig,Et,Esig);
    end
    FWHMt = fwhm(magsq(Et),tau);
    Ew = fftc(Et);    
    FWHMw = fwhm(magsq(Ew),WaveLen);    
    fprintf('Iteration = %d\t\tFWHM_t = %.2ffs\t\tFWHM_w = %.2fnm\t\tTime_sec = %.2f\n',...
        j,FWHMt,FWHMw,toc) 
    
    if(PlotIter < MaxIter)
        close(hhh);
    end
    
end

