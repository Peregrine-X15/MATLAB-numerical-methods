% clear all
% clc
% close all

%for courser part 81->41, and 41->21
N=41;
M=41;
N_c=21;
M_c=21;
Res_Cond=1e-6;
del_x=1/(N-1);
del_y=1/(M-1);
del_xc=1/(N_c-1);
del_yc=1/(M_c-1);

x=linspace(0,1,N);
y=linspace(0,1,M);

count_GMG_small=0;
R2_GMG=1;
phi_g=zeros(N,M);
phi_c=zeros(N_c,M_c);
phi_corr=zeros(N,M);
ResO=zeros(N,M);
Res_C=zeros(N_c,M_c);
phi_g(N,:)=100.*(1-y)+500*exp(-50.*y.^2);
phi_g(1,:)=500*exp(-50*(1+y.^2));
phi_g(:,1)=100*x+500*exp(-50*(1-x).^2);
phi_g(:,M)=500*exp(-50*((1-x).^2+1));
%making fine mesh
for i=2:N-1
    for j=2:M-1
        S_phi(i,j)=50000.*exp(-50.*((1-x(i)).^2+y(j).^2)).*(100.*((1-x(i)).^2+y(j).^2)-2);
    end
end

while R2_GMG>Res_Cond && count_GMG_small<20000
for i=2:N-1
    for j=2:M-1
              %phi equation from presentation, Gauss Seidel update for GMG
         phi_g(i,j)=(phi_g(i+1,j)/del_x^2+phi_g(i-1,j)/del_x^2 ...
             +phi_g(i,j-1)/del_y^2 +phi_g(i,j+1)/del_y^2-S_phi(i,j))/(2/del_y^2+2/del_x^2);
    end
end

   count_GMG_small=count_GMG_small+1;
R2_GMG=0;
for i=2:N-1
    for j=2:M-1
%   Residual calculation
% RHS and LHS, consult page 137 of pdf, 127 of book 
ResO(i,j)=S_phi(i,j)-phi_g(i+1,j)/del_x^2-phi_g(i-1,j)/del_x^2 ...
             -phi_g(i,j-1)/del_y^2 -phi_g(i,j+1)/del_y^2+phi_g(i,j)*(2/del_y^2+2/del_x^2);
R2_GMG=R2_GMG+ResO(i,j)*ResO(i,j);
    end
end
R2_GMG=sqrt(R2_GMG);
R2_vectGMG_small(count_GMG_small)=R2_GMG;

%transfering redisdual into course mesh interpolation
for I=1:N_c
    for J=1:M_c
        i = 2*I-1;
        j=2*J-1;
        Res_C(I,J)=ResO(i,j);
    end
end

for I=2:N_c-1
    for J=2:M_c-1
        phi_c(I,J)=Res_C(I,J)/(-2/del_xc^2-2/del_yc^2);
    end
end
%Fine mesh correction by interpolation
for I=1:N_c
    for J=1:M_c
        i = 2*I-1;
        j= 2*J-1;
        %direct transfer
        phi_corr(i,j)=phi_c(I,J);
    end
end
for I=1:N_c-1
    for J=1:M_c-1
        i = 2*I-1;
        j= 2*J-1;
        %Two point
        %page 204 in book
        phi_corr(i+1,j)=(phi_c(I,J)+phi_c(I+1,J))/2;
        phi_corr(i,j+1)=(phi_c(I,J)+phi_c(I,J+1))/2;
        %Four point
        %page 205 in book
        phi_corr(i+1,j+1)=(phi_c(I,J)+phi_c(I+1,J)+phi_c(I,J+1)+phi_c(I+1,J+1))/4;
    end
end
phi_g=phi_g+phi_corr;
end  

phi_g=phi_g';
for i=1:M/2
t_fact=phi_g(M+1-i,:);
phi_g(M+1-i,:)=phi_g(i,:);
phi_g(i,:)=t_fact;
end
%phi_star=phi_star';
[X, Y] = meshgrid(x,y(end:-1:1));
%calculate analytical here.
phi_ana=500*exp(-50*((1-X).^2+Y.^2))+100*X.*(1-Y);
%check equations for errors.
%transpose, by adding ' to phi_star and switching rows
% figure(1)
% contour(X,Y,phi_g,'ShowText','on','LevelStep',50)
% xlabel('x')
% ylabel('y')
% 
% %R2, error with analyitical contour
% figure (2)
% semilogy(1:count_GMG_small,R2_vectGMG_small(1:count_GMG_small))
% % filename='Residuals_41.xlsx';
% % xlswrite([filename],[(1:count)' R2_vect(1:count)'],'Sheet1');
% % xlswrite([filename],[(1:count_GS)' R2_vectGS(1:count_GS)'],'Sheet2');
% xlabel('iterations')
% ylabel('Residual')
% legend('GMG')
% hold off
% 
% figure (3)
% contour(X,Y,phi_ana-phi_g,'ShowText','on','LevelStep',0.1)
% xlabel('x')
% ylabel('y')