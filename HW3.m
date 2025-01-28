clear all
close all
clc
%% Jacobi method
N=41;
M=41;
Res_Cond=1e-6;
del_x=1/(N-1);
del_y=1/(M-1);

x=linspace(0,1,N);
y=linspace(0,1,M);

count=0;
R2=1;
phi_star=zeros(N,M);
phi_g=zeros(N,M);


while R2>Res_Cond && count<20000
for i=2:N-1
    for j=2:M-1
        %k=(j-1)*N+i;
          %phi_star(i,j)=zeros(N,M);

    %BCs (check if eqns are correct)
        phi_star(N,j)=100.*(1-y(j))+500*exp(-50.*y(j).^2);
        phi_star(1,j)=500*exp(-50*(1+y(j)^2));
        phi_star(i,1)=100*x(i)+500*exp(-50*(1-x(i))^2);
        phi_star(i,M)=500*exp(-50*((1-x(i)).^2+1));

        phi_g(N,j)=100.*(1-y(j))+500*exp(-50.*y(j).^2);
        phi_g(1,j)=500*exp(-50*(1+y(j)^2));
        phi_g(i,1)=100*x(i)+500*exp(-50*(1-x(i))^2);
        phi_g(i,M)=500*exp(-50*((1-x(i)).^2+1));

%S_phi
S_phi(i,j)=50000.*exp(-50.*((1-x(i)).^2+y(j).^2)).*(100.*((1-x(i)).^2+y(j).^2)-2);

        %phi equation from presentation
         phi_g(i,j)=(phi_star(i+1,j)/del_x^2+phi_star(i-1,j)/del_x^2 ...
             +phi_star(i,j-1)/del_y^2 +phi_star(i,j+1)/del_y^2-S_phi(i,j))/(2/del_y^2+2/del_x^2);
       
    end
end
count=count+1;
phi_star=phi_g;
R2=0;
for i=2:N-1
    for j=2:M-1
%   Residual calculation
% RHS and LHS, consult page 137 of pdf, 127 of book 
ResO=S_phi(i,j)-phi_star(i+1,j)/del_x^2-phi_star(i-1,j)/del_x^2 ...
             -phi_star(i,j-1)/del_y^2 -phi_star(i,j+1)/del_y^2+phi_star(i,j)*(2/del_y^2+2/del_x^2);
R2=R2+ResO*ResO;
    end
end
R2=sqrt(R2);
R2_vect(count)=R2;
end  

%transposing
phi_star=phi_star';
for i=1:M/2
t_fact=phi_star(M+1-i,:);
phi_star(M+1-i,:)=phi_star(i,:);
phi_star(i,:)=t_fact;
end
%phi_star=phi_star';
[X, Y] = meshgrid(x,y(end:-1:1));
%calculate analytical here.
phi_ana=500*exp(-50*((1-X).^2+Y.^2))+100*X.*(1-Y);
%check equations for errors.
%transpose, by adding ' to phi_star and switching rows
figure(1)
contour(X,Y,phi_star,'ShowText','on')
xlabel('x')
ylabel('y')

%R2, error with analyitical contour
figure (2)
semilogy(1:count,R2_vect(1:count))
xlabel('iterations')
ylabel('Residual')

figure (3)
contour(X,Y,phi_ana-phi_star,'ShowText','on')
xlabel('x')
ylabel('y')

%% Gauss-Seidel
N=41;
M=41;
Res_Cond=1e-6;
del_x=1/(N-1);
del_y=1/(M-1);

x=linspace(0,1,N);
y=linspace(0,1,M);

count_GS=0;
R2_GS=1;
phi_star=zeros(N,M);
phi_g=zeros(N,M);


while R2_GS>Res_Cond && count_GS<20000
for i=2:N-1
    for j=2:M-1
        %k=(j-1)*N+i;
          %phi_star(i,j)=zeros(N,M);

    %BCs (check if eqns are correct)
        phi_star(N,j)=100.*(1-y(j))+500*exp(-50.*y(j).^2);
        phi_star(1,j)=500*exp(-50*(1+y(j)^2));
        phi_star(i,1)=100*x(i)+500*exp(-50*(1-x(i))^2);
        phi_star(i,M)=500*exp(-50*((1-x(i)).^2+1));

        phi_g(N,j)=100.*(1-y(j))+500*exp(-50.*y(j).^2);
        phi_g(1,j)=500*exp(-50*(1+y(j)^2));
        phi_g(i,1)=100*x(i)+500*exp(-50*(1-x(i))^2);
        phi_g(i,M)=500*exp(-50*((1-x(i)).^2+1));

%S_phi
S_phi(i,j)=50000.*exp(-50.*((1-x(i)).^2+y(j).^2)).*(100.*((1-x(i)).^2+y(j).^2)-2);

        %phi equation from presentation
         phi_g(i,j)=(phi_star(i+1,j)/del_x^2+phi_g(i-1,j)/del_x^2 ...
             +phi_g(i,j-1)/del_y^2 +phi_star(i,j+1)/del_y^2-S_phi(i,j))/(2/del_y^2+2/del_x^2);
       
    end
end
count_GS=count_GS+1;
phi_star=phi_g;
R2_GS=0;
for i=2:N-1
    for j=2:M-1
%   Residual calculation
% RHS and LHS, consult page 137 of pdf, 127 of book 
ResO=S_phi(i,j)-phi_star(i+1,j)/del_x^2-phi_star(i-1,j)/del_x^2 ...
             -phi_star(i,j-1)/del_y^2 -phi_star(i,j+1)/del_y^2+phi_star(i,j)*(2/del_y^2+2/del_x^2);
R2_GS=R2_GS+ResO*ResO;
    end
end
R2_GS=sqrt(R2_GS);
R2_vectGS(count_GS)=R2_GS;
end  

%transposing
phi_star=phi_star';
for i=1:M/2
t_fact=phi_star(M+1-i,:);
phi_star(M+1-i,:)=phi_star(i,:);
phi_star(i,:)=t_fact;
end
%phi_star=phi_star';
[X, Y] = meshgrid(x,y(end:-1:1));
%calculate analytical here.
phi_ana=500*exp(-50*((1-X).^2+Y.^2))+100*X.*(1-Y);
%check equations for errors.
%transpose, by adding ' to phi_star and switching rows
figure(4)
contour(X,Y,phi_star,'ShowText','on')
xlabel('x')
ylabel('y')

%R2, error with analyitical contour
figure (5)
semilogy(1:count,R2_vect(1:count))
hold on
semilogy(1:count_GS,R2_vectGS(1:count_GS))
% filename='Residuals_41.xlsx';
% xlswrite([filename],[(1:count)' R2_vect(1:count)'],'Sheet1');
% xlswrite([filename],[(1:count_GS)' R2_vectGS(1:count_GS)'],'Sheet2');
xlabel('iterations')
ylabel('Residual')
legend('Jacobi','Gauss-Seidel')
hold off

figure (6)
contour(X,Y,phi_ana-phi_star,'ShowText','on','LevelStep',0.1)
xlabel('x')
ylabel('y')