clear all
close all
clc
%for 40 x 40
tic
N=40;
M=40;
Res_Cond=1e-6;
del_x=1/(N);
del_y=1/(M);

x=del_x/2:del_x:1-del_x/2;
y=del_y/2:del_y:1-del_y/2;
% x(1,40)=1;
% y(1,40)=1;
[X, Y] = meshgrid(x,y(end:-1:1));
phi_ana=500*exp(-50*((1-X).^2+Y.^2))+100*X.*(1-Y);
for i=1:N
    for j=1:M
        S_phi(i,j)=50000.*exp(-50.*((1-x(i)).^2+y(j).^2)).*(100.*((1-x(i)).^2+y(j).^2)-2);
    end
end


    %BCs 
        phi_E=100.*(1-y)+500*exp(-50.*y.^2);
        phi_W=500*exp(-50*(1+y.^2));
        phi_S=100.*x+500*exp(-50.*(1-x).^2);
        phi_N=500*exp(-50*((1-x).^2+1));
%phi_e=phi(i+1,j)
%phi_W=phi(i-1,j)
%phi_n=phi(i,j+1)
%phi_s=phi(i,j-1)
phi=zeros(N,M);
res=zeros(N,M);
R2=1;
count=0;
while R2>Res_Cond && count<10000
for i=1:N
    for j=1:M
        if i~=1 && i~=N && j~=1 && j~=M
        %Interior
        phi(i,j)=(-S_phi(i,j)+(1/del_x^2)*phi(i+1,j)+(1/del_x^2)*phi(i-1,j)+(1/del_y^2)*phi(i,j+1)+1/(del_y^2)*phi(i,j-1))/(2/del_x^2+2/del_y^2);
        end
%west side
if i==1  && j~=1 && j~=M
phi(i,j)=(S_phi(i,j)*del_x*del_y-4*del_y/(3*del_x)*phi(i+1,j)-8*del_y/(3*del_x)*phi_W(j)-(del_x/del_y)*phi(i,j+1)-del_x/del_y*phi(i,j-1))/(-4*del_y/del_x-2*del_x/del_y);
end
%East side
if i==N  && j~=1 && j~=M
phi(i,j)=(S_phi(i,j)*del_y*del_x-4*del_y/(3*del_x)*phi(i-1,j)-8*del_y/(3*del_x)*phi_E(j)-del_x/del_y*phi(i,j+1)-del_x/del_y*phi(i,j-1))/(-4*del_y/del_x-2*del_x/del_y);
end
%North
if j==M && i~=1 && i~=N
    phi(i,j)=(S_phi(i,j)*del_y*del_x-4*del_x/(3*del_y)*phi(i,j-1)-8*del_x/(3*del_y)*phi_N(i)-del_y/del_x*phi(i+1,j)-del_y/del_x*phi(i-1,j))/(-4*del_x/del_y-2*del_y/del_x);
end
%South
if j==1 && i~=1 && i~=N
phi(i,j)=(S_phi(i,j)*del_y*del_x-del_y/del_x*phi(i+1,j)-4*del_x/(3*del_y)*phi(i,j+1)-8*del_x/(3*del_y)*phi_S(i)-del_y/del_x*phi(i-1,j))/(-2*del_x/del_y-4*del_y/del_x);
end       
%Northeast
if j==M && i==N
phi(i,j)=(S_phi(i,j)*del_x*del_y-4*del_y/(3*del_x)*phi(i-1,j)-8*del_y/(3*del_x)*phi_E(j)-4*del_x/(3*del_y)*phi(i,j-1)-8*del_x/(3*del_y)*phi_N(i))/(-4*del_y/del_x-4*del_x/del_y);
end
%Southeast
%
if i==N && j==1
phi(i,j)=(S_phi(i,j)*del_x*del_y-8*del_y/(3*del_x)*phi_E(j)-8*del_x/(3*del_y)*phi_S(i)-4*del_y/(3*del_x)*phi(i-1,j)-4*del_x/(3*del_y)*phi(i,j+1))/(-4*del_y/del_x-4*del_x/del_y);
end
%southwest
if i==1 && j==1
% phi(i,j)=(S_phi(i,j)*del_x*del_y-5*del_y/(3*del_x)*phi(i+1,j)-8*del_x/(3*del_y)*phi_S(i)-del_y/(3*del_x)*phi_W(j)-5*del_x/(3*del_y)*phi(i,j+1))/(2*del_y/del_x-4*del_x/del_y);
phi(i,j)=(S_phi(i,j)*del_x*del_y-4*del_y/(3*del_x)*phi(i+1,j)-8*del_x/(3*del_y)*phi_S(i)-8*del_y/(3*del_x)*phi_W(j)-4*del_x/(3*del_y)*phi(i,j+1))/(-4*del_y/del_x-4*del_x/del_y);
end
%Northwest
if i==1 && j==M
phi(i,j)=(S_phi(i,j)*del_y*del_x-4*del_y/(3*del_x)*phi(i+1,j)-8*del_y/(3*del_x)*phi_W(j)-4*del_x/(3*del_y)*phi(i,j-1)-8*del_x/(3*del_y)*phi_N(i))/(-4*del_y/del_x-4*del_x/del_y);
end
    end
end
count=count+1;
for i=1:N
    for j=1:M
% Residuals
if i~=1 && i~=N && j~=1 && j~=M
    %Interior
    res(i,j)=phi(i,j)-(-S_phi(i,j)+(1/del_x^2)*phi(i+1,j)+(1/del_x^2)*phi(i-1,j)+(1/del_y^2)*phi(i,j+1)+1/(del_y^2)*phi(i,j-1))/(2/del_x^2+2/del_y^2);
end
%west side
if i==1  && j~=1 && j~=M
res(i,j)=phi(i,j)-(S_phi(i,j)*del_x*del_y-4*del_y/(3*del_x)*phi(i+1,j)-8*del_y/(3*del_x)*phi_W(j)-(del_x/del_y)*phi(i,j+1)-del_x/del_y*phi(i,j-1))/(-4*del_y/del_x-2*del_x/del_y);
end
%East side
if i==N  && j~=1 && j~=M
res(i,j)=phi(i,j)-(S_phi(i,j)*del_y*del_x-4*del_y/(3*del_x)*phi(i-1,j)-8*del_y/(3*del_x)*phi_E(j)-del_x/del_y*phi(i,j+1)-del_x/del_y*phi(i,j-1))/(-4*del_y/del_x-2*del_x/del_y);
end
%North
if j==M && i~=1 && i~=N
    res(i,j)=phi(i,j)-(S_phi(i,j)*del_y*del_x-4*del_x/(3*del_y)*phi(i,j-1)-8*del_x/(3*del_y)*phi_N(i)-del_y/del_x*phi(i+1,j)-del_y/del_x*phi(i-1,j))/(-4*del_x/del_y-2*del_y/del_x);
end
%South
if j==1 && i~=1 && i~=N
res(i,j)=phi(i,j)-(S_phi(i,j)*del_y*del_x-del_y/del_x*phi(i+1,j)-4*del_x/(3*del_y)*phi(i,j+1)-8*del_x/(3*del_y)*phi_S(i)-del_y/del_x*phi(i-1,j))/(-2*del_x/del_y-4*del_y/del_x);
end       
%Northeast
if j==M && i==N
res(i,j)=phi(i,j)-(S_phi(i,j)*del_x*del_y-4*del_y/(3*del_x)*phi(i-1,j)-8*del_y/(3*del_x)*phi_E(j)-4*del_x/(3*del_y)*phi(i,j-1)-8*del_x/(3*del_y)*phi_N(i))/(-4*del_y/del_x-4*del_x/del_y);
end
%Southeast
if i==N && j==1
res(i,j)=phi(i,j)-(S_phi(i,j)*del_x*del_y-8*del_y/(3*del_x)*phi_E(j)-8*del_x/(3*del_y)*phi_S(i)-4*del_y/(3*del_x)*phi(i-1,j)-4*del_x/(3*del_y)*phi(i,j+1))/(-4*del_y/del_x-4*del_x/del_y);
end
%southwest
if i==1 && j==1
% res(i,j)=phi(i,j)-(S_phi(i,j)*del_x*del_y-5*del_y/(3*del_x)*phi(i+1,j)-8*del_x/(3*del_y)*phi_S(i)-del_y/(3*del_x)*phi_W(j)-5*del_x/(3*del_y)*phi(i,j+1))/(2*del_y/del_x-4*del_x/del_y);
res(i,j)=phi(i,j)-(S_phi(i,j)*del_x*del_y-4*del_y/(3*del_x)*phi(i+1,j)-8*del_x/(3*del_y)*phi_S(i)-8*del_y/(3*del_x)*phi_W(j)-4*del_x/(3*del_y)*phi(i,j+1))/(-4*del_y/del_x-4*del_x/del_y);
end
%Northwest
if i==1 && j==M
res(i,j)=phi(i,j)-(S_phi(i,j)*del_y*del_x-4*del_y/(3*del_x)*phi(i+1,j)-8*del_y/(3*del_x)*phi_W(j)-4*del_x/(3*del_y)*phi(i,j-1)-8*del_x/(3*del_y)*phi_N(i))/(-4*del_y/del_x-4*del_x/del_y);
end
    end
end
% res_squared=res.^2;
% R2=sum(res_squared,'all');
% R2=sqrt(R2);

%     end
% end
res_squared=res.^2;
R2=sum(res_squared,'all');
R2=sqrt(R2);
end
% Residual calculation
% Copy paste above code but make adjustment as below
%res(i,j) = phi(i,j)-(-S_phi(i,j)+(1/del_x^2)*phi(i+1,j)+(1/del_x^2)*phi(i-1,j)+(1/del_y^2)*phi(i,j+1)+1/(del_y^2)*phi(i,j-1))/(2/del_x^2+2/del_y^2);
% Rest of copied code
% res_squared = res.^2;
% R2 = sum(res_squared,'all');
% R2=sqrt(R2);

%fluxes
%treat same as FDM
%phi_e=phi(i+1,j)
%phi_W=phi(i-1,j)
%phi_n=phi(i,j+1)
%phi_s=phi(i,j-1)
for i=1:N
    %added delx and dely to all Qs
    for j=1:M
        if j==M && i~=1 && i~=N
            %north
       %delphi_dely_N=(9*phi(i,j)-phi_S(i)-8*phi_N(i))/(-3*del_y);
       delphi_dely_N=(9*phi(i,j)-phi(i,j-1)-8*phi_N(i))/(-3*del_y);
       QN=delphi_dely_N*del_x*del_y;
        end
        if j==1 && i~=1 && i~=N
            %south
       %delphi_dely_S=(9*phi(i,j)-phi_N(i)-8*phi_S(i))/(3*del_y);
       delphi_dely_S=(9*phi(i,j)-phi(i,j+1)-8*phi_S(i))/(3*del_y);
       QS=delphi_dely_S*del_x*del_y; 
        end
        if i==N  && j~=1 && j~=M
            %east
       %delphi_delx_E=(9*phi(i,j)-phi_W(j)-8*phi_E(j))/(-3*del_x);
       delphi_delx_E=(9*phi(i,j)-phi(i-1,j)-8*phi_E(j))/(-3*del_x);
       QE=delphi_delx_E*del_y*del_x;
        end
        if i==1  && j~=1 && j~=M
            %west
       %delphi_delx_W=(9*phi(i,j)-phi_E(j)-8*phi_W(j))/(3*del_x);
       delphi_delx_W=(9*phi(i,j)-phi(i+1,j)-8*phi_W(j))/(3*del_x);
       QW=delphi_delx_W*del_y*del_x;
        end
    end
end
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
            P(i,j)=S_phi(i,j)*del_x*del_y;
        end
        %east side
        if i==N
            P(i,j)=S_phi(i,j)*del_x*del_y;
        end
          %South Side
        if j==1
            P(i,j)=S_phi(i,j)*del_x*del_y;
        end
        %North Side
        if j==M
            P(i,j)=S_phi(i,j)*del_x*del_y;
        end
            %Corners
        %south west corner
        if i==1 && j==1
            P(i,j)=S_phi(i,j)*del_x*del_y;    
        end
        %south east corner
        if i==1 && j==M
            P(i,j)=S_phi(i,j)*del_x*del_y; 
        end
        %North east corner
        if i==N && j==M
            P(i,j)=S_phi(i,j)*del_x*del_y; 
        end
        %North West corner
        if i==N && j==1
            P(i,j)=S_phi(i,j)*del_x*del_y; 
        end
    end
end
P=sum(P,'all');

I=QE-QW+QN-QS-P;
%compute residuals, based on same eqns as coded in.
phi=phi';
for i=1:M/2
t_fact=phi(M+1-i,:);
phi(M+1-i,:)=phi(i,:);
phi(i,:)=t_fact;
end
phi=phi';
figure 
contour(X,Y,phi_ana-phi,'ShowText','on','LevelStep',5)
xlabel('x')
ylabel('y')
figure
contour(X,Y,phi,'ShowText','on')
xlabel('X')
ylabel('Y')
toc