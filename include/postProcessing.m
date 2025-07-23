% Script postProcessing.m processes the raw signal from the specified data
% file (in settings) operating on blocks of 37 seconds of data.
%
% First it runs acquisition code identifying the satellites in the file,
% then the code and carrier for each of the satellites are tracked, storing
% the 1msec accumulations.  After processing all satellites in the 37 sec
% data block, then postNavigation is called. It calculates pseudoranges
% and attempts a position solutions. At the end plots are made for that
% block of data.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis, Dennis M. Akos
% Some ideas by Dennis M. Akos
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%                         THE SCRIPT "RECIPE"
%
% The purpose of this script is to combine all parts of the software
% receiver.
%
% 1.1) Open the data file for the processing and seek to desired point.
%
% 2.1) Acquire satellites
%
% 3.1) Initialize channels (preRun.m).
% 3.2) Pass the channel structure and the file identifier to the tracking
% function. It will read and process the data. The tracking results are
% stored in the trackResults structure. The results can be accessed this
% way (the results are stored each millisecond):
% trackResults(channelNumber).XXX(fromMillisecond : toMillisecond), where
% XXX is a field name of the result (e.g. I_P, codePhase etc.)
%
% 4) Pass tracking results to the navigation solution function. It will
% decode navigation messages, find satellite positions, measure
% pseudoranges and find receiver position.
%
% 5) Plot the results.


%% Initialization =========================================================
disp ('   Starting processing...');

[fid, message] = fopen(settings.fileName, 'rb');

%Initialize the multiplier to adjust for the data type
if (settings.fileType==1)
    dataAdaptCoeff=1;
else
    dataAdaptCoeff=2;
end

%If success, then process the data
if (fid > 0)

    % Move the starting point of processing. Can be used to start the
    % signal processing at any point in the data record (e.g. good for long
    % records or for signal processing in blocks).
    fseek(fid, dataAdaptCoeff*settings.skipNumberOfSamples, 'bof');

    %% Acquisition ============================================================

    % Do acquisition if it is not disabled in settings or if the variable
    % acqResults does not exist.
    if ((settings.skipAcquisition == 0))
        % Find number of samples per spreading code
        samplesPerCode = round(settings.samplingFreq / ...
            (settings.codeFreqBasis / settings.codeLength));
        %--- Do the acquisition -------------------------------------------
        disp ('   Acquiring satellites...');
        % Read the required amount of data depending on the data file type
        % and the number of code period of coherent and non-coherent 
        % integration and invoke the acquisition function
        data = fread(fid, dataAdaptCoeff*(settings.acquisition.cohCodePeriods* ...
            settings.acquisition.nonCohSums+1)*samplesPerCode*6, settings.dataType)';
        if (dataAdaptCoeff==2)
            data1=data(1:2:end);
            data2=data(2:2:end);
            data=data1 + 1i .* data2;
        end
        acqResults = acquisition(data, settings);
       % Plot the acquisition results
        plotAcquisition(acqResults);
    else
        disp('   skipAcquistion==1');
        acqResults = load([settings.dir '\acqResults.mat']);
        plotAcquisition(acqResults);
    end

    %% Initialize channels and prepare for the run ============================

    % Start further processing only if a GNSS signal was acquired (the
    % field FREQUENCY will be set to 0 for all not acquired signals)
    if (any(acqResults.peakMetric>settings.acqThreshold))
        channel = preRun(acqResults, settings);
        showChannelStatus(channel, settings);
    else
        % No satellites to track, exit
        disp('   No GNSS signals detected, signal processing finished.');
        trackResults = [];
        return;
    end

    %% Track the signal =======================================================
    if ~settings.INSAT
        startTime = now;
        disp (['   Tracking started at ', datestr(startTime)]);
        [trackResults, channel] = tracking(fid, channel, settings);
        % Close the data file
        fclose(fid);
        disp(['   Tracking is over (elapsed time ',datestr(now - startTime, 13), ')'] )
        %Compute the PRM C/No and add it to the track results
        trackResults = calculateCNoPRM(trackResults,settings);
        %Compute the MOM C/No and add it to the track results
        trackResults = calculateCNoMOM(trackResults,settings);
        % Auto save the acquisition & tracking results to a file to allow
        % running the positioning solution afterwards.
        disp('   Saving Acq & Tracking results to file "trackingResults.mat"')
        save([settings.dir '\acqResults.mat'], 'acqResults');
        save([settings.dir '\trkResults.mat'], 'trackResults', 'settings', 'acqResults', 'channel');
    else
        disp('   skip scalar tracking, load exisiting results')
        trackResults = load([settings.dir '\trkResults.mat']).trkResults;
    end

    %% Calculate navigation solutions =========================================
    
    if ~settings.INSAT
        disp('   Calculating navigation solutions...');
        settings = initSettings();
        trackResults = load([settings.dir '\trkResults.mat']).trackResults;
        [navSolutions, eph, svTimeTable,activeChnList] = postNavigationNomial(trackResults, settings);
        save([settings.dir '\navSolutions.mat'],'navSolutions','eph','svTimeTable','activeChnList')
    else
        disp('   skip postNavigation, load existing navSolution')
        %input the path and file name of the existing navSolution saved
        %from scalar tracking and calculation
        settings = initSettings();
        [fid, message] = fopen(settings.fileName, 'rb');
        load([settings.dir '\navSolutions.mat']);
        disp('   start INSAT')
        %start vector tracking
        [Rx,trackResults] = INSAT(fid, channel,trackResults,navSolutions,eph,activeChnList,svTimeTable, settings);
        imuInfo = load('E:\TEXBAT\postResult\cleanDynamic\navResults.mat').navResults;
        save([settings.dir '\INSAT.mat'],'Rx','trackResults')
        navSolutions = load([settings.dir '\navSolutions.mat']).navSolutions;
        load([settings.dir '\INSAT.mat'])
        plotINSAT(Rx,imuInfo,navSolutions,trackResults)
        fclose all;
        return;     
    end
disp('   Processing is complete for this data block');

else
    % Error while opening the data file.
    error('Unable to read file %s: %s.', settings.fileName, message);
end % if (fid > 0)

