function [trackResults]= calculateCNoPRM(trackResults,settings)

%Calculate CNo using the PRM Method
%
%[trackResults]= calculateCNoPRM(trackResults,settings)
%
%   Inputs:
%       trackResults    - Results from Tracking
%       settings        - Settings Structure
%   Outputs:
%       trackResults    - trackResuls appended with C/No data
%
%--------------------------------------------------------------------------
% Copyright (C) D.M.Akos
% Written by Sirish Jetti
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

%For each channel
for channelNr = 1:size(trackResults,2)

    %Find the bit start position
    navBits=trackResults(channelNr).I_P;
    navBits(navBits>0)=1;
    navBits(navBits<=0)=-1;
    navBitStartPosIndex = find(diff(navBits),10);
    navBitStartPos = 0;
    for i=2:numel(navBitStartPosIndex)

        navBitStartPos=navBitStartPosIndex(i);
        if ((navBitStartPosIndex(i)-navBitStartPosIndex(i-1))...
                ==20/(1000*settings.CNo.accTime))
            break;
        end
    end


    navBitStartPos = mod(navBitStartPos,20/(1000*settings.CNo.accTime))+1;

    %Compute the C/No(PRM)
    trackResults(channelNr).CNo.PRMValue=CNoPRM(...
        trackResults(channelNr).I_P(1,navBitStartPos:end),...
        trackResults(channelNr).Q_P(1,navBitStartPos:end),...
        settings.CNo.PRM_K,settings.CNo.PRM_M,settings.CNo.accTime);

    %Compute the time index for calculated C/No
    trackResults(channelNr).CNo.PRMIndex=...
        navBitStartPos:settings.CNo.PRM_K:...
        navBitStartPos+(numel(trackResults(channelNr).CNo.PRMValue)-1)*settings.CNo.PRM_K;
end