function plotINSAT(Rx,imuInfo,navSolutions,trackResults,settings)
% Achieve the initialization of the receiver and INSAT
%
%   Inputs:
%       Rx              - initialization of the receiver
%       imuInfo         - authenic information
%       navSolusions    - navigation solutions from normal receiver
%       trackResults    - INSAT results

    %% plot trajectory
    plotLen = round((length(Rx.est_lat*180/pi)+1000)/500);
    figure,geoplot(imuInfo.INSlat(1:plotLen*50),imuInfo.INSlon(1:plotLen*50),'LineWidth',3)
    hold on
    geoplot(navSolutions.latitude,navSolutions.longitude,'LineWidth',3)
    geoplot(Rx.est_lat*180/pi,Rx.est_lon*180/pi,'-.','LineWidth',3)
    legend('Authentic','Spoofing','INSAT','Location','northwest')
    
    %% plot Multi-correlator 
    if settings.multiCorrOn
        figure
        iChannel = 1;
        trackRsultFlipMultiCorr = zeros(size(trackResults(iChannel).multiCorr));
        FlipOrNot = max(abs(trackResults(iChannel).multiCorr))+min(trackResults(iChannel).multiCorr)>0;
        for ilen = 1:length(trackResults(iChannel).multiCorr)  
            trackRsultFlipMultiCorr(:,ilen) = trackResults(iChannel).multiCorr(:,ilen)*FlipOrNot(ilen);
        end
        trackRsultFlipMultiCorr(trackRsultFlipMultiCorr<-1e5) = 0;
        [X,Y] = meshgrid(0.2:0.2:300,settings.multiCorrDistribution);
        surf(X,Y,trackRsultFlipMultiCorr(:,1:200:3e5))
        ylabel('Code Phase(Chips)'),xlabel('Times(s)'),zlabel('Correlation Value'),view([108 36]) %view([112 45])    
    end    
end