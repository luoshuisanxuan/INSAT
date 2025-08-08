function [Rx,INSATsetting] = initINSATsetting(settings,channel,trackRes,svTimeTable,activeChnList,navSolutions)
% Achieve the initialization of the receiver and INSAT
%
%   Inputs:
%       fid             - file identifier of the signal record for iloopCnt+1
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRum.iChannel from
%                       acquisition results).
%       trackRes        -tracking results using scalar loop to initialize
%       navSolusions    -navigation solutions from scalar loop to
%                       initialize
%       eph             - ephemerides
%       activeChnList   - a list of active satellites in the dataset
%       svTimeTable     - satellite time to find transmit time of a sample 
%       settings        - receiver settings.
%   Outputs:
%       Rx              - initialization of the receiver
%       INSATsetting    - initialization of INSAT

    %% base parameter
    INSATsetting.kpt = 1e-3;           % kalman pprocess time /s
    INSATsetting.kmt = 1e-3;           % kalman measurement time /s
    INSATsetting.pdi = 1e-3;           % integration time /s
    INSATsetting.r2d = 180/pi;         % rad to deg
    INSATsetting.D2r = pi/180;         % deg to rad
    INSATsetting.StartTime = 1000;     % start time for INSAT /ms
    INSATsetting.tracklength = 400000; % total track length /ms
    %adaptive filtering window for the Kalman filter measurements
    INSATsetting.cnt=1;
    INSATsetting.lastn=50;
    %% INS Measurements  
    texbatCleanDynamic = load('texbatCleanDynamic.mat');
    INSATsetting.GPSTime = texbatCleanDynamic.GPSTime';        
    INSATsetting.INSlat = texbatCleanDynamic.INSlat';
    INSATsetting.INSlon = texbatCleanDynamic.INSlon';
    INSATsetting.INShei = texbatCleanDynamic.INShei';
    INSATsetting.INSroll = texbatCleanDynamic.INSroll;
    INSATsetting.INSpitch = texbatCleanDynamic.INSpitch;
    INSATsetting.INShead = texbatCleanDynamic.INShead;
    INSATsetting.INSve = texbatCleanDynamic.INSve';
    INSATsetting.INSvn = texbatCleanDynamic.INSvn';
    INSATsetting.INSvu = texbatCleanDynamic.INSvu';
    INSATsetting.INSaccy = texbatCleanDynamic.INSaccy;
    INSATsetting.INSaccx = texbatCleanDynamic.INSaccx;
    INSATsetting.INSaccz = texbatCleanDynamic.INSaccz;
    INSATsetting.INSgyroy = texbatCleanDynamic.INSgyroy;
    INSATsetting.INSgyrox = texbatCleanDynamic.INSgyrox;
    INSATsetting.INSgyroz = texbatCleanDynamic.INSgyroz;
    INSATsetting.INSaccby = texbatCleanDynamic.INSaccby;
    INSATsetting.INSaccbx = texbatCleanDynamic.INSaccbx;
    INSATsetting.INSaccbz = texbatCleanDynamic.INSaccbz;
    INSATsetting.INSgyrody = texbatCleanDynamic.INSgyrody;
    INSATsetting.INSgyrodx = texbatCleanDynamic.INSgyrodx;
    INSATsetting.INSgyrodz = texbatCleanDynamic.INSgyrodz;
    
    %% KF parament
    INSATsetting.stateno = 17; %number of states
    INSATsetting.Qw = diag([diag(1e0*eye(3))',diag(1e-3*eye(3))',1*diag(1e-2*eye(3))',1*diag(1e-8*eye(3))',1*diag(1e-8*eye(3))',1e-6,1e-1]);
    
    %measurement noise var-covariance matrix
    INSATsetting.R(1:settings.numberOfChannels,1:settings.numberOfChannels) = 1500*eye(settings.numberOfChannels);
    INSATsetting.R(settings.numberOfChannels+1:2*settings.numberOfChannels,settings.numberOfChannels+1:2*settings. ...
        numberOfChannels) = 9e2*eye(settings.numberOfChannels);
    % initial estimation error var-covairance matrix
    INSATsetting.P0 = diag([1e0,1e0,1e0,1e-1,1e-1,1e-1,1*diag(1e-10*eye(3))',1*diag(1e-10*eye(3))',1*diag(1e-10*eye(3))',1,1e-8]);
    
    %initialize measurement matrix
    INSATsetting.H = zeros(2*settings.numberOfChannels,INSATsetting.stateno);
    
    %initialize measurement vector
    INSATsetting.Z = zeros(2*settings.numberOfChannels,1);
    % states of Kalman filter initialization
    INSATsetting.X_est = zeros(INSATsetting.stateno,INSATsetting.tracklength);
    INSATsetting.deltaX = zeros(1,INSATsetting.tracklength);
    INSATsetting.alphak = zeros(1,INSATsetting.tracklength);
    INSATsetting.X0 = zeros(INSATsetting.stateno,1);
    
    %% Receiver Initial
    npts = INSATsetting.tracklength/10+1; %number of points in IMU measurement dataset
    %find the true position for the INSATsetting.StartTime sample point
    ind1 = find(INSATsetting.GPSTime>=navSolutions.rxTime(INSATsetting.StartTime*settings.navSolRate/1000),1);
    lat0 = INSATsetting.INSlat(ind1)*INSATsetting.D2r;
    lon0 = INSATsetting.INSlon(ind1)*INSATsetting.D2r;
    hei0 = INSATsetting.INShei(ind1);
    [pos0(1,1),pos0(1,2),pos0(1,3)] = geo2cart([lat0*INSATsetting.r2d,0,0],[lon0*INSATsetting.r2d,0,0], hei0, 5);
    Rx.pos_kf = pos0;
    %find attitude for the INSATsetting.StartTime sample point
    ind0 = find(INSATsetting.GPSTime>=navSolutions.rxTime(INSATsetting.StartTime*settings.navSolRate/1000),1);   
    phi = INSATsetting.INSroll(ind0)/INSATsetting.r2d;
    theta = INSATsetting.INSpitch(ind0)/INSATsetting.r2d;
    psi = INSATsetting.INShead(ind0)/INSATsetting.r2d;
    %direction cosine matrix
    DCMnb = eul2dcm([phi theta psi]);
    %initialize the output to be saved
    Rx.est_roll_KF = zeros(1,INSATsetting.tracklength+1);
    Rx.est_pitch_KF = zeros(1,INSATsetting.tracklength+1);
    Rx.est_yaw_KF = zeros(1,INSATsetting.tracklength+1);
    Rx.est_roll_KF(1) = phi;
    Rx.est_pitch_KF(1) = theta;
    Rx.est_yaw_KF(1) = psi;
    Rx.ve = INSATsetting.INSve(ind0);
    Rx.vn = INSATsetting.INSvn(ind0);
    Rx.vu = INSATsetting.INSvu(ind0);
    %initialize intermediate variables for INS update
    [tlat,tlon,thei] = cart2geo(pos0(1,1),pos0(1,2),pos0(1,3),5);
    orginllh = [tlat*INSATsetting.D2r,tlon*INSATsetting.D2r,thei];
    Rx.est_lat = zeros(1,INSATsetting.tracklength+1);
    Rx.est_lon = zeros(1,INSATsetting.tracklength+1);
    Rx.est_height = zeros(1,INSATsetting.tracklength+1);
    Rx.est_lat(1) = orginllh(1);
    Rx.est_lat(2) = orginllh(1);
    Rx.est_lon(1) = orginllh(2);
    Rx.est_height(1) = orginllh(3);
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
    slat = sin(Rx.est_lat(1));   clat = cos(Rx.est_lat(1));
    slon = sin(Rx.est_lon(1));   clon = cos(Rx.est_lon(1));
    est_DCMel=[-slon -slat*clon clat*clon
            clon -slat*slon clat*slon
            0    clat       slat]';
    Rx.est_DCMel_KF = est_DCMel;
    Rx.omega_e = 7.292115e-5;
    Rx.omega_ie_E = [0 0 Rx.omega_e]';
    Rx.C = [0 1 0; 1 0 0; 0 0 -1];  
    %transform IMU raw measurements to delta velocity and delta theta in 1ms
    ind0 = find(INSATsetting.GPSTime>=navSolutions.rxTime(INSATsetting.StartTime*settings.navSolRate/1000),1);
    Rx.raw_dv = 0.001*[INSATsetting.INSaccy(ind0:ind0+npts-1)-0*INSATsetting.INSaccby(ind1:ind1+npts-1),INSATsetting.INSaccx(ind0:ind0+npts-1)...
        -0*INSATsetting.INSaccbx(ind1:ind1+npts-1),-INSATsetting.INSaccz(ind0:ind0+npts-1)+0*INSATsetting.INSaccbz(ind1:ind1+npts-1)];
    Rx.raw_dtheta = 0.001*INSATsetting.D2r*[INSATsetting.INSgyroy(ind0:ind0+npts-1)-0*INSATsetting.INSgyrody(ind1:ind1+npts-1),INSATsetting.INSgyrox(ind0:ind0+npts-1)...
        -0*INSATsetting.INSgyrodx(ind1:ind1+npts-1),-INSATsetting.INSgyroz(ind0:ind0+npts-1)+0*INSATsetting.INSgyrodz(ind1:ind1+npts-1)];
    
    %% Receiver INS Aided Tracking Initial
    Rx.mat1 = zeros(settings.numberOfChannels,INSATsetting.lastn);
    Rx.mat2 = zeros(settings.numberOfChannels,INSATsetting.lastn);
    %initialize intermediate variables
    Rx.dSv = zeros(1,settings.numberOfChannels);
    Rx.dPlos = zeros(1,settings.numberOfChannels);
    Rx.Vs = zeros(1,settings.numberOfChannels);
    Rx.dVlos = zeros(1,settings.numberOfChannels);
    Rx.carrError = zeros(1,settings.numberOfChannels);
    Rx.carrErrorold = zeros(1,settings.numberOfChannels);
    Rx.codeError = zeros(1,settings.numberOfChannels);
    Rx.codeErrorold = zeros(1,settings.numberOfChannels);
    Rx.carrFreq = zeros(1,settings.numberOfChannels);
    Rx.codeFreq = zeros(1,settings.numberOfChannels);
    initsample = zeros(1,settings.numberOfChannels);
    initsampleforcode = zeros(1,settings.numberOfChannels);
    Rx.codePhase = zeros(1,settings.numberOfChannels);
    Rx.codePhaseStep = zeros(1,settings.numberOfChannels);
    Rx.carrFreqBasis = zeros(1,settings.numberOfChannels);
    Rx.remCarrPhase = zeros(1,settings.numberOfChannels);
    Rx.remCodePhase = zeros(1,settings.numberOfChannels);
    Rx.oldCodeNco = zeros(1,settings.numberOfChannels);
    Rx.oldCodeError = zeros(1,settings.numberOfChannels);
    Rx.oldCarrNco = zeros(1,settings.numberOfChannels);
    Rx.oldCarrError = zeros(1,settings.numberOfChannels);
    Rx.carrNco = zeros(1,settings.numberOfChannels);
    Rx.vsmCnt = zeros(1,settings.numberOfChannels);
    Rx.CNo = 0;
    Rx.satPosenu = zeros(3,settings.numberOfChannels);
    Rx.satPosenu0 = zeros(3,settings.numberOfChannels);
    Rx.satVelenu = zeros(3,settings.numberOfChannels);
    Rx.blksize = zeros(1,settings.numberOfChannels);
    %initialize code, frequency, transmit time
    for channelNr = 1:settings.numberOfChannels
            Rx.carrFreq(1,channelNr) = trackRes(1,activeChnList(channelNr)).carrFreq(INSATsetting.StartTime);
            Rx.carrFreqBasis(1,channelNr) = channel(channelNr).acquiredFreq;
            Rx.codeFreq(1,channelNr) = trackRes(1,activeChnList(channelNr)).codeFreq(INSATsetting.StartTime);
            initsample(1,channelNr) = ceil(trackRes(1,activeChnList(channelNr)).absoluteSample(INSATsetting.StartTime));
            initsampleforcode(1,channelNr) = ceil(trackRes(1,activeChnList(channelNr)).absoluteSample(INSATsetting.StartTime-1));
            Rx.codePhase(1,channelNr) = (initsampleforcode(1,channelNr)-trackRes(1,activeChnList(channelNr)).absoluteSample ...
                (INSATsetting.StartTime-1))/settings.samplingFreq*Rx.codeFreq(1,channelNr);
            Rx.codePhaseStep(1,channelNr) = Rx.codeFreq(1,channelNr) / settings.samplingFreq;
            tTime = findTransTime(initsample(channelNr),activeChnList,svTimeTable,trackRes);
            Rx.transmitTime(activeChnList(channelNr)) = tTime(activeChnList(channelNr));
            Rx.remCarrPhase(1,channelNr) = trackRes(1,activeChnList(channelNr)).remCarrPhase(INSATsetting.StartTime);
            Rx.remCodePhase(1,channelNr) = trackRes(1,activeChnList(channelNr)).remCodePhase(INSATsetting.StartTime);
            Rx.oldCodeNco(1,channelNr) = trackRes(1,activeChnList(channelNr)).dllDiscrFilt(INSATsetting.StartTime-1);
            Rx.oldCodeError(1,channelNr) = trackRes(1,activeChnList(channelNr)).dllDiscr(INSATsetting.StartTime-1);
            Rx.oldCarrNco(1,channelNr) = trackRes(1,activeChnList(channelNr)).pllDiscrFilt(INSATsetting.StartTime-1);
            Rx.oldCarrError(1,channelNr) = trackRes(1,activeChnList(channelNr)).pllDiscr(INSATsetting.StartTime-1);
            %C/No computation
            Rx.vsmCnt(channelNr)  = 0;
            caCode0 = generateCAcode(trackRes(1,activeChnList(channelNr)).PRN);
            if settings.multiCorrOn
                Rx.caCode(channelNr,:) =[caCode0(1020) caCode0(1021) caCode0(1022) caCode0(1023) caCode0 caCode0(1) caCode0(2) caCode0(3) caCode0(4)];
            else
                Rx.caCode(channelNr,:) =[caCode0(1023) caCode0 caCode0(1)];
            end
            Rx.blksize(1,channelNr) = ceil((settings.codeLength-Rx.remCodePhase(1,channelNr)) / Rx.codePhaseStep(1,channelNr));
    end % for channelNr
    Rx.transmitTime0 = Rx.transmitTime;
    Rx.blksize0 = INSATsetting.pdi*settings.samplingFreq*ones(1,settings.numberOfChannels);
    Rx.samplepos = initsample;
    mininit = min(initsample);
    Rx.minpos = mininit;
    Rx.IP1 = zeros(1,settings.numberOfChannels);
    Rx.QP1 = zeros(1,settings.numberOfChannels);
end