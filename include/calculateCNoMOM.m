function [trackResults]= calculateCNoMOM(trackResults,settings)
%Calculate CNo using the Moments Method
%
%[trackResults]= calculateCNoMOM(trackResults,settings)
%
%   Inputs:
%       trackResults    - Results from Tracking
%       settings        - Settings Structure
%   Outputs:
%       trackResults    - trackResuls appended with C/No data
%
%--------------------------------------------------------------------------
% Copyright (C) D.M.Akos
% Written by Xiaofan Li
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
%
%USA.
%--------------------------------------------------------------------------


%For each channel
for channelNr = 1:size(trackResults,2)
    % Calculate the total number of CNo estimation
    numCNo=floor(settings.msToProcess/settings.ms/settings.CNo.MOMinterval);
    % Compute CNo at each measuring point
    for momCnt=1:numCNo
        %Compute the time index for calculated C/No
        loopCnt=momCnt*settings.CNo.MOMinterval;
        CNoValue=CNoMOM(trackResults(channelNr).I_P(loopCnt-settings.CNo.MOMinterval+1:loopCnt),...
                        trackResults(channelNr).Q_P(loopCnt-settings.CNo.MOMinterval+1:loopCnt),settings.CNo.accTime);
        trackResults(channelNr).CNo.MOMValue(momCnt)=CNoValue;            
        trackResults(channelNr).CNo.MOMIndex(momCnt)=loopCnt;
    end
end;

