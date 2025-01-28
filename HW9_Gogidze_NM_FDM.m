clear all
clc
close all
tic
HW3;
%for 41x41
%uses same difference plot as in HW3
for i=1:N
    for j=1:M
        S_phi(i,j)=50000.*exp(-50.*((1-x(i)).^2+y(j).^2)).*(100.*((1-x(i)).^2+y(j).^2)-2);
    end
end
% QE=zeros(1,M);
% QS=zeros(1,N);
%Boundary Conditions
        for j=1:M
        %WEST BC
        delPhi_delx_W=(4*phi_star(1+1,j)-phi_star(1+2,j)-3*phi_star(1,j))/(2*del_x);
         %west
        QW(1,j)=delPhi_delx_W*del_y;
        end
        for j=1:M
        %EAST BC
        % delPhi_delx_E(j)=-(4*phi_star(N,j)-phi_star(N-2,j)-3*phi_star(N,j))/(2*del_x);
        delPhi_delx_E(j)=-(4*phi_star(j,N)-phi_star(j,N-2)-3*phi_star(j,N))/(2*del_x);
        QE(1,j)=delPhi_delx_E(j)*del_y;
        end
        for i=1:N
        %NORTH BC
        delPhi_dely_N=-(4*phi_star(i,end-1)-phi_star(i,end-2)-3*phi_star(i,end))/(2*del_y);
        QN(1,i)=delPhi_dely_N*del_x;
        end
        for i=1:N
        %SOUTH BC
        % delPhi_dely_S(i)=(4*phi_star(i,1)-phi_star(i,1+2)-3*phi_star(i,1))/(2*del_y);
         delPhi_dely_S=(4*phi_star(i,1+1)-phi_star(i,1+2)-3*phi_star(i,1))/(2*del_y);
        QS(1,i)=delPhi_dely_S*del_x; 
        end
       
%first corners
QE(1)=QE(1,1)*0.5;
QW(1)=QW(1,1)*0.5;
QS(1)=QS(1,1)*0.5;
QN(1)=QN(1,1)*0.5;
%more corners
QE(end)=QE(1,end)*0.5;
QW(end)=QW(1,end)*0.5;
QS(end)=QS(1,end)*0.5;
QN(end)=QN(1,end)*0.5;

        QE=sum(QE);
        QW=sum(QW);
        QN=sum(QN);
        QS=sum(QS);

P=zeros(41,41);
for i=1:N
    for j=1:M
        P(i,j)=S_phi(i,j)*del_x*del_y;
        %Sides
       %west side 
        if i==1
            P(i,j)=S_phi(i,j)*del_x*del_y*0.5;
        end
        %east side
        if i==N
            P(i,j)=S_phi(i,j)*del_x*del_y*0.5;
        end
          %South Side
        if j==1
            P(i,j)=S_phi(i,j)*del_x*del_y*0.5;
        end
        %North Side
        if j==M
            P(i,j)=S_phi(i,j)*del_x*del_y*0.5;
        end
            %Corners
        %south west corner
        if i==1 && j==1
            P(i,j)=S_phi(i,j)*del_x*del_y*0.25;    
        end
        %south east corner
        if i==1 && j==M
            P(i,j)=S_phi(i,j)*del_x*del_y*0.25; 
        end
        %North east corner
        if i==N && j==M
            P(i,j)=S_phi(i,j)*del_x*del_y*0.25; 
        end
        %North West corner
        if i==N && j==1
            P(i,j)=S_phi(i,j)*del_x*del_y*0.25; 
        end
    end
end
P=sum(P,'all');
%look over Qs
%all qs are summed not minus
I=QE-QW+QN-QS-P;
toc
%half for all corners then sum all together.

%plots are from HW3