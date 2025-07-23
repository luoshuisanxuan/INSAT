function [F,radiusa] = Fupdate(Rx,INSATsetting,iloopCnt)
% Update function for the error propagation coefficient matrix
%
%   Inputs:
%       Rx              - initialization of the receiver
%       INSATsetting    - initialization of INSAT
%       iloopCnt        - epoch
%   Outputs:
%       F              - initialization of the receiver
%       radiusa        - earth radius
    
    %% parameter calculation
    radiusa = localradius(Rx.est_lat(iloopCnt+1));
    slat = sin(Rx.est_lat(iloopCnt+1));
    clat = cos(Rx.est_lat(iloopCnt+1));
    Rx.ve = Rx.vel_l(iloopCnt+1,1);
    Rx.vn = Rx.vel_l(iloopCnt+1,2);
    Rx.vu = Rx.vel_l(iloopCnt+1,3);
    rh = radiusa+Rx.est_height(iloopCnt+1);
    %% error propagation coefficient matrix
    F11 = zeros(3,3);
    F11(1,2) = Rx.ve*slat/rh/clat/clat;
    F11(1,3) = Rx.ve/rh/rh/clat;
    F11(2,3) = Rx.vn/rh/rh;
    F12 = eye(3);
    F13 = zeros(3,3);   
    F21=zeros(3,3);
    F21(1,2) = -Rx.ve*(2*Rx.omega_e*clat+Rx.ve/rh/clat/clat);
    F21(1,3) = (Rx.ve*Rx.ve*slat/clat+Rx.vn*Rx.vu)/rh/rh;
    F21(2,2) = 2*Rx.omega_e*(Rx.vn*clat+Rx.vu*slat)+Rx.vn*Rx.ve/rh/clat/clat;
    F21(2,3) = -Rx.ve/rh/rh*(Rx.vn*slat/clat-Rx.vu);
    F21(3,2) = -2*Rx.omega_e*Rx.ve*slat;
    F21(3,3) = (Rx.ve*Rx.ve+Rx.vn*Rx.vn)/rh/rh;    
    F22 = zeros(3,3);
    F22(1,1) = (Rx.vn*slat/clat-Rx.vu)/rh;
    F22(1,2) = 2*Rx.omega_e*slat+Rx.ve*slat/clat/rh;
    F22(1,3) = -(2*Rx.omega_e*clat+Rx.ve/rh);
    F22(2,1) = -2*(Rx.omega_e*slat+Rx.ve*slat/clat/rh);
    F22(2,2) = -Rx.vu/rh;
    F22(2,3) = -Rx.vn/rh;
    F22(3,1) = 2*(Rx.omega_e*clat+Rx.ve/rh);
    F22(3,2) = 2*Rx.vn/rh;    
    F23 = antisymm(Rx.accel_L);  
    F31 = zeros(3,3);
    F32 = zeros(3,3);
    F33 = (-1)*antisymm(Rx.omega_il_L);      
    F = zeros(17,17);
    F(1:9,1:9) = [F11 F12 F13; F21 F22 F23; F31 F32 F33]; 
    F(10:12,10:12) = zeros(3,3);
    F(4:6,10:12) = zeros(3,3);   
    F(13:15,13:15) = zeros(3,3);
    F(7:9,13:15) = Rx.C*Rx.est_DCMbn;
    F(16,17) = 1;
    % discretize F matrix, assuming all states are irrelavent to
    % themselves at last epoch, except for clock bias and drift
    F = F*INSATsetting.kmt+eye(INSATsetting.stateno);
    F(1,1) = 0;
    F(2,2) = 0;
    F(3,3) = 0;
    F(4,4) = 0;
    F(5,5) = 0;
    F(6,6) = 0;
    F(7,7) = 0;
    F(8,8) = 0;
    F(9,9) = 0;
    F(16,16) = 1;
    F(17,17) = 0;
end