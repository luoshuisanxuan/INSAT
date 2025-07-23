function plotINSAT(Rx,imuInfo,navSolutions,trackResults)
% Achieve the initialization of the receiver and INSAT
%
%   Inputs:
%       Rx              - initialization of the receiver
%       imuInfo         - authenic information
%       navSolusions    - navigation solutions from normal receiver
%       trackResults    - INSAT results

    %% plot trajectory
    plotLen = round((length(Rx.est_lat*180/pi)+1000)/500);
    figure,geoplot(imuInfo.latitude(1:plotLen),imuInfo.longitude(1:plotLen),'LineWidth',3)
    hold on
    geoplot(navSolutions.latitude,navSolutions.longitude,'LineWidth',3)
    geoplot(Rx.est_lat*180/pi,Rx.est_lon*180/pi,'-.','LineWidth',3)
    legend('Authentic','Spoofing','INSAT','Location','northwest')
    
    %% plot Multi-correlator 
    figure
    MultiCorrDistribution=-3:0.1:3;
    trkResults = trackResults;
    iChannel = 1;
    trackRsultFlipMultiCorr = zeros(size(trkResults(iChannel).multiCorr));
    FlipOrNot = max(abs(trkResults(iChannel).multiCorr))+min(trkResults(iChannel).multiCorr)>0;
    for ilen = 1:length(trkResults(iChannel).multiCorr)  
        trackRsultFlipMultiCorr(:,ilen) = trkResults(iChannel).multiCorr(:,ilen)*FlipOrNot(ilen);
    end
    trackRsultFlipMultiCorr(trackRsultFlipMultiCorr<-1e5) = 0;
    [X,Y] = meshgrid(0.2:0.2:300,MultiCorrDistribution);
    surf(X,Y,trackRsultFlipMultiCorr(:,1:200:3e5))
    ylabel('Code Phase(Chips)'),xlabel('Times(s)'),zlabel('Correlation Value'),view([108 36]) %view([112 45])    
end