function [CNo]= CNoMOM(I,Q,T)
%Calculate CNo using the Moments method
%
%[CNo]= CNoMOM(I,Q,T)
%
%   Inputs:
%       I           - Prompt In Phase values of the signal from Tracking
%       Q           - Prompt Quadrature Phase values of the signal from Tracking
%       T          - Accumulation interval in Tracking (in sec)
%   Outputs:
%       CNo         - Estimated C/No for the given values of I and Q
%
%
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
%USA.
%--------------------------------------------------------------------------

% Compute the normalized noise equivelent bandwidth
Beq=1/T;

% Second moment 
m2=mean(abs(I+j*Q).^2);
% Fourth moment
m4=mean(abs(I+j*Q).^4);
% Signal level
Ps=sqrt(2*m2^2-m4);
% Noise level
Pn=m2-Ps;
if Pn>0
    %Calculate C/No
    CNo=10*log10(Ps/Pn*Beq);
else
    CNo=0;
end

