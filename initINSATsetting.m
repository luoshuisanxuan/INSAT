function [Rx,INSATsetting] = initINSATsetting(settings,channel,trackRes,svTimeTable,activeChnList,navSolutions)
    %% base parameter
    INSATsetting.kpt = 1e-3;%kalman filter process update time
    INSATsetting.kmt = 1e-3;%kalman filter measurement update time/s
    INSATsetting.pdi = 1e-3;%integration time /s
    INSATsetting.r2d = 180/pi;
    INSATsetting.D2r = pi/180;
    INSATsetting.StartTime = 1000;%start time in dataset /ms
    INSATsetting.begintime = INSATsetting.StartTime/1000;% start time for plots in the end /s
    INSATsetting.tracklengthall = 36000;  %TOTAL TRACKING LENGTH in ms
    INSATsetting.tracklength = 36000; % vector tracking length FOR ONE SEGMENT ms
    %adaptive filtering window for the Kalman filter measurements
    INSATsetting.cnt=1;
    INSATsetting.lastn=50;
    %% INS Measurements
    imuInfo = load('F:\TEXBAT\postResult\cleanDynamic\navSolutions.mat').navSolutions;
    wgs84 = wgs84Ellipsoid('meter');
    Ve = zeros(1,length(imuInfo.X));
    Vn = zeros(1,length(imuInfo.X));
    Vu = zeros(1,length(imuInfo.X));
    for i = 1:length(imuInfo.X)
        [Ve1,Vn1,Vu1] = ecef2enu(imuInfo.Vx(i)*0.1+imuInfo.X(i),imuInfo.Vy(i)*0.1+imuInfo.Y(i), ...
            imuInfo.Vz(i)*0.1+imuInfo.Z(i),imuInfo.latitude(i),imuInfo.longitude(i),imuInfo.height(i),wgs84);
        [Ve2,Vn2,Vu2] = ecef2enu(imuInfo.X(i),imuInfo.Y(i),imuInfo.Z(i), ...
            imuInfo.latitude(i),imuInfo.longitude(i),imuInfo.height(i),wgs84);
        Ve(i) = (Ve1-Ve2)/0.1;Vn(i) = (Vn1-Vn2)/0.1;Vu(i) = (Vu1-Vu2)/0.1;
    end
    lenNav = 4098;
    INSATsetting.GPSTime = (navSolutions.rxTime(1):0.01:navSolutions.rxTime(end))';        
    INSATsetting.INSlat = interp1(1:lenNav,imuInfo.latitude,1:0.1:lenNav)';
    INSATsetting.INSlon = interp1(1:lenNav,imuInfo.longitude,1:0.1:lenNav)';
    INSATsetting.INShei = interp1(1:lenNav,imuInfo.height,1:0.1:lenNav)';
    INSATsetting.INSroll = zeros(size(INSATsetting.INShei));
    INSATsetting.INSpitch = zeros(size(INSATsetting.INShei));
    INSATsetting.INShead = zeros(size(INSATsetting.INShei));
    INSATsetting.INSve = interp1(1:lenNav,Ve,1:0.1:lenNav)';
    INSATsetting.INSvn = interp1(1:lenNav,Vn,1:0.1:lenNav)';
    INSATsetting.INSvu = interp1(1:lenNav,Vu,1:0.1:lenNav)';
    INSATsetting.INSaccy = zeros(size(INSATsetting.INShei));
    INSATsetting.INSaccx = zeros(size(INSATsetting.INShei));
    INSATsetting.INSaccz = zeros(size(INSATsetting.INShei));
    INSATsetting.INSgyroy = zeros(size(INSATsetting.INShei));
    INSATsetting.INSgyrox = zeros(size(INSATsetting.INShei));
    INSATsetting.INSgyroz = zeros(size(INSATsetting.INShei));
    INSATsetting.INSaccby = zeros(size(INSATsetting.INShei));
    INSATsetting.INSaccbx = zeros(size(INSATsetting.INShei));
    INSATsetting.INSaccbz = zeros(size(INSATsetting.INShei));
    INSATsetting.INSgyrody = zeros(size(INSATsetting.INShei));
    INSATsetting.INSgyrodx = zeros(size(INSATsetting.INShei));
    INSATsetting.INSgyrodz = zeros(size(INSATsetting.INShei));
    
    %% KF parament
    INSATsetting.stateno=17; %number of states
    %process noise var-covariance matrix
    INSATsetting.Qw=diag([diag(1e0*eye(3))',diag(1e-3*eye(3))',1*diag(1e-2*eye(3))',1*diag(1e-8*eye(3))',1*diag(1e-8*eye(3))',1e-6,1e-1]);
    
    %measurement noise var-covariance matrix
    INSATsetting.R(1:settings.numberOfChannels,1:settings.numberOfChannels)=1500*eye(settings.numberOfChannels);
    INSATsetting.R(settings.numberOfChannels+1:2*settings.numberOfChannels,settings.numberOfChannels+1:2*settings. ...
        numberOfChannels)=9e2*eye(settings.numberOfChannels);
    % initial estimation error var-covairance matrix
    INSATsetting.P0=diag([1e0,1e0,1e0,1e-1,1e-1,1e-1,1*diag(1e-10*eye(3))',1*diag(1e-10*eye(3))',1*diag(1e-10*eye(3))',1,1e-8]);
    
    %initialize measurement matrix
    INSATsetting.H=zeros(2*settings.numberOfChannels,INSATsetting.stateno);
    
    %initialize measurement vector
    INSATsetting.Z=zeros(2*settings.numberOfChannels,1);
    % states of Kalman filter initialization
    INSATsetting.X_est=zeros(INSATsetting.stateno,INSATsetting.tracklength);
    INSATsetting.deltaX = zeros(1,INSATsetting.tracklength);
    INSATsetting.alphak =  zeros(1,INSATsetting.tracklength);
    INSATsetting.X0=zeros(INSATsetting.stateno,1);
    
    %% Receiver Initial
    npts = INSATsetting.tracklengthall/10+1;%number of points in IMU measurement dataset
    %find the true position for the INSATsetting.StartTime sample point
    ind1=find(INSATsetting.GPSTime>=navSolutions.rxTime(INSATsetting.StartTime*settings.navSolRate/1000),1);
    lat0=INSATsetting.INSlat(ind1)*INSATsetting.D2r;
    lon0=INSATsetting.INSlon(ind1)*INSATsetting.D2r;
    hei0=INSATsetting.INShei(ind1);
    [pos0(1,1),pos0(1,2),pos0(1,3)]=geo2cart([lat0*INSATsetting.r2d,0,0],[lon0*INSATsetting.r2d,0,0], hei0, 5);
    Rx.pos_kf=pos0;
    %find attitude for the INSATsetting.StartTime sample point
    ind0=find(INSATsetting.GPSTime>=navSolutions.rxTime(INSATsetting.StartTime*settings.navSolRate/1000),1);   
    phi=INSATsetting.INSroll(ind0)/INSATsetting.r2d;
    theta=INSATsetting.INSpitch(ind0)/INSATsetting.r2d;
    psi=INSATsetting.INShead(ind0)/INSATsetting.r2d;
    %direction cosine matrix
    DCMnb=eul2dcm([phi theta psi]);
    %initialize the output to be saved
    Rx.est_roll_KF=zeros(1,INSATsetting.tracklength+1);
    Rx.est_pitch_KF=zeros(1,INSATsetting.tracklength+1);
    Rx.est_yaw_KF=zeros(1,INSATsetting.tracklength+1);
    Rx.est_roll_KF(1) = phi;
    Rx.est_pitch_KF(1) = theta;
    Rx.est_yaw_KF(1) = psi;
    Rx.ve=INSATsetting.INSve(ind0);
    Rx.vn=INSATsetting.INSvn(ind0);
    Rx.vu=INSATsetting.INSvu(ind0);
    %initialize intermediate variables for INS update
    [tlat,tlon,thei]=cart2geo(pos0(1,1),pos0(1,2),pos0(1,3),5);
    orginllh=[tlat*INSATsetting.D2r,tlon*INSATsetting.D2r,thei];
    Rx.est_lat=zeros(1,INSATsetting.tracklength+1);
    Rx.est_lon=zeros(1,INSATsetting.tracklength+1);
    Rx.est_height=zeros(1,INSATsetting.tracklength+1);
    Rx.est_lat(1)=orginllh(1);
    Rx.est_lat(2)=orginllh(1);
    Rx.est_lon(1)=orginllh(2);
    Rx.est_height(1)=orginllh(3);
    height = orginllh(3); 
    
    Rx.heightold = height;
    Rx.veold = Rx.ve;
    Rx.vnold = Rx.vn;
    Rx.vuold = Rx.vu;
    Rx.vel_l(1,:) = [Rx.veold Rx.vnold Rx.vu];
    Rx.velenu=Rx.vel_l(1,:);
    Rx.velold = [Rx.ve, Rx.vn, Rx.vu];
    Rx.latold = orginllh(1);
    Rx.est_DCMbn = DCMnb';
    Rx.est_DCMbn_KF = Rx.est_DCMbn;
    slat=sin(Rx.est_lat(1));   clat=cos(Rx.est_lat(1));
    slon=sin(Rx.est_lon(1));   clon=cos(Rx.est_lon(1));
    est_DCMel=[-slon -slat*clon clat*clon
            clon -slat*slon clat*slon
            0    clat       slat]';
    Rx.est_DCMel_KF=est_DCMel;
    Rx.omega_e=7.292115e-5;
    Rx.omega_ie_E=[0 0 Rx.omega_e]';
    Rx.C = [0 1 0; 1 0 0; 0 0 -1];    % conversion between NED and ENU
    %transform IMU raw measurements to delta velocity and delta theta in 1ms
    ind0=find(INSATsetting.GPSTime>=navSolutions.rxTime(INSATsetting.StartTime*settings.navSolRate/1000),1);
    Rx.raw_dv=0.001*[INSATsetting.INSaccy(ind0:ind0+npts-1)-0*INSATsetting.INSaccby(ind1:ind1+npts-1),INSATsetting.INSaccx(ind0:ind0+npts-1)...
        -0*INSATsetting.INSaccbx(ind1:ind1+npts-1),-INSATsetting.INSaccz(ind0:ind0+npts-1)+0*INSATsetting.INSaccbz(ind1:ind1+npts-1)];
    Rx.raw_dtheta=0.001*INSATsetting.D2r*[INSATsetting.INSgyroy(ind0:ind0+npts-1)-0*INSATsetting.INSgyrody(ind1:ind1+npts-1),INSATsetting.INSgyrox(ind0:ind0+npts-1)...
        -0*INSATsetting.INSgyrodx(ind1:ind1+npts-1),-INSATsetting.INSgyroz(ind0:ind0+npts-1)+0*INSATsetting.INSgyrodz(ind1:ind1+npts-1)];
    
    %% Receiver INS Aided Tracking Initial
    Rx.mat1=zeros(settings.numberOfChannels,INSATsetting.lastn);
    Rx.mat2=zeros(settings.numberOfChannels,INSATsetting.lastn);
    %initialize intermediate variables
    Rx.dSv=zeros(1,settings.numberOfChannels);
    Rx.dPlos=zeros(1,settings.numberOfChannels);
    Rx.Vs=zeros(1,settings.numberOfChannels);
    Rx.dVlos=zeros(1,settings.numberOfChannels);
    Rx.carrError=zeros(1,settings.numberOfChannels);
    Rx.carrErrorold=zeros(1,settings.numberOfChannels);
    Rx.codeError=zeros(1,settings.numberOfChannels);
    Rx.codeErrorold=zeros(1,settings.numberOfChannels);
    Rx.carrFreq=zeros(1,settings.numberOfChannels);
    Rx.codeFreq=zeros(1,settings.numberOfChannels);
    initsample=zeros(1,settings.numberOfChannels);
    initsampleforcode=zeros(1,settings.numberOfChannels);
    Rx.codePhase=zeros(1,settings.numberOfChannels);
    Rx.codePhaseStep=zeros(1,settings.numberOfChannels);
    Rx.carrFreqBasis=zeros(1,settings.numberOfChannels);
    Rx.remCarrPhase=zeros(1,settings.numberOfChannels);
    Rx.remCodePhase=zeros(1,settings.numberOfChannels);
    Rx.oldCodeNco=zeros(1,settings.numberOfChannels);
    Rx.oldCodeError=zeros(1,settings.numberOfChannels);
    Rx.oldCarrNco=zeros(1,settings.numberOfChannels);
    Rx.oldCarrError=zeros(1,settings.numberOfChannels);
    Rx.carrNco=zeros(1,settings.numberOfChannels);
    Rx.vsmCnt=zeros(1,settings.numberOfChannels);
    Rx.CNo=0;
    Rx.satPosenu = zeros(3,settings.numberOfChannels);
    Rx.satPosenu0 = zeros(3,settings.numberOfChannels);
    Rx.satVelenu = zeros(3,settings.numberOfChannels);
    Rx.blksize=zeros(1,settings.numberOfChannels);
    %initialize code, frequency, transmit time
    for channelNr = 1:settings.numberOfChannels%settings.numberOfChannels    
        % Only process if PRN is non zero (acquisition was successful)
            % Rx.trackResults(activeChnList(channelNr)).PRN     = trackRes(1,activeChnList(channelNr)).PRN;
            Rx.carrFreq(1,channelNr)=trackRes(1,activeChnList(channelNr)).carrFreq(INSATsetting.StartTime);
            Rx.carrFreqBasis(1,channelNr) = channel(channelNr).acquiredFreq;
            Rx.codeFreq(1,channelNr)=trackRes(1,activeChnList(channelNr)).codeFreq(INSATsetting.StartTime);
            initsample(1,channelNr)=ceil(trackRes(1,activeChnList(channelNr)).absoluteSample(INSATsetting.StartTime));
            initsampleforcode(1,channelNr)=ceil(trackRes(1,activeChnList(channelNr)).absoluteSample(INSATsetting.StartTime-1));
            Rx.codePhase(1,channelNr)=(initsampleforcode(1,channelNr)-trackRes(1,activeChnList(channelNr)).absoluteSample ...
                (INSATsetting.StartTime-1))/settings.samplingFreq*Rx.codeFreq(1,channelNr);
            Rx.codePhaseStep(1,channelNr) = Rx.codeFreq(1,channelNr) / settings.samplingFreq;
            tTime=findTransTime(initsample(channelNr),activeChnList,svTimeTable,trackRes);
            Rx.transmitTime(activeChnList(channelNr))=tTime(activeChnList(channelNr));
            Rx.remCarrPhase(1,channelNr)=trackRes(1,activeChnList(channelNr)).remCarrPhase(INSATsetting.StartTime);
            Rx.remCodePhase(1,channelNr)=trackRes(1,activeChnList(channelNr)).remCodePhase(INSATsetting.StartTime);
            Rx.oldCodeNco(1,channelNr)=trackRes(1,activeChnList(channelNr)).dllDiscrFilt(INSATsetting.StartTime-1);
            Rx.oldCodeError(1,channelNr)=trackRes(1,activeChnList(channelNr)).dllDiscr(INSATsetting.StartTime-1);
            Rx.oldCarrNco(1,channelNr)=trackRes(1,activeChnList(channelNr)).pllDiscrFilt(INSATsetting.StartTime-1);
            Rx.oldCarrError(1,channelNr)=trackRes(1,activeChnList(channelNr)).pllDiscr(INSATsetting.StartTime-1);
            %C/No computation
            Rx.vsmCnt(channelNr)  = 0;
            % Get a vector with the C/A code sampled 1x/chip
            caCode0 = generateCAcode(trackRes(1,activeChnList(channelNr)).PRN);
            Rx.caCode(channelNr,:) =[caCode0(1023) caCode0 caCode0(1)];
            Rx.blksize(1,channelNr)=ceil((settings.codeLength-Rx.remCodePhase(1,channelNr)) / Rx.codePhaseStep(1,channelNr));
    end % for channelNr
    Rx.transmitTime0=Rx.transmitTime;
    Rx.blksize0=INSATsetting.pdi*settings.samplingFreq*ones(1,settings.numberOfChannels);
    Rx.samplepos=initsample;
    mininit=min(initsample);
    Rx.minpos=mininit;
    Rx.IP1=zeros(1,settings.numberOfChannels);
    Rx.QP1=zeros(1,settings.numberOfChannels);
end