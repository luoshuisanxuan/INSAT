function plotINSAT(Rx,imuInfo,navSolutions)
    figure,geoplot(imuInfo.latitude,imuInfo.longitude,'LineWidth',4)
    hold on
    geoplot(navSolutions.latitude,navSolutions.longitude)
    geoplot(Rx.est_lat*180/pi,Rx.est_lon*180/pi,'LineWidth',4)
end
%% Multi-correlator 
% figure
% MultiCorrDistribution=-3:0.1:3;
% trkResults = trackResults;
% iChannel = 1;
% trackRsultFlipMultiCorr = zeros(size(trkResults(iChannel).multiCorr));
% FlipOrNot = max(abs(trkResults(iChannel).multiCorr))+min(trkResults(iChannel).multiCorr)>0;
% for ilen = 1:length(trkResults(iChannel).multiCorr)  
%     trackRsultFlipMultiCorr(:,ilen) = trkResults(iChannel).multiCorr(:,ilen)*FlipOrNot(ilen);
% end
% trackRsultFlipMultiCorr(trackRsultFlipMultiCorr<-1e5) = 0;
% surf(trackRsultFlipMultiCorr(:,1:200:3e5))