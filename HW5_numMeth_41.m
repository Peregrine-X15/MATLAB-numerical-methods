clear all
close all
clc
%% Method of Steepest Descent
% N=[41 81];
% M=[41 81];
N=41;
M=41;
neq=N*M;
phi_g=zeros(1,neq);
Res_cond=10^-6;
%intial variables

del_x=1/(N-1);
del_y=1/(M-1);
Res_Cond=1e-6;
x=linspace(0,1,N);
y=linspace(0,1,M);
%BC's
phi_R=100.*(1-y)+500*exp(-50.*y.^2);
phi_L=500*exp(-50*(1+y.^2));
phi_B=100*x+500*exp(-50*(1-x).^2);
phi_T=500*exp(-50*((1-x).^2+1));
%diagonals from HW4
EE=(2/(del_x^2)+2/(del_y^2))*ones(1,neq);
BB=[zeros(1,N),(-1/del_y^2)*ones(1,neq-N)];
DD=[0,-1/(del_x^2)*ones(1,neq-1)];
FF=[(-1/(del_x^2))*ones(1,neq-1),0];
HH=[(-1/(del_y^2))*ones(1,neq-N),zeros(1,N)];

    for i=1:N
       for j=1:M
       k=(j-1)*N+i;
       QQ(k)=-50000.*exp(-50.*((1-x(i)).^2+y(j).^2)).*(100.*((1-x(i)).^2+y(j).^2)-2);
  
       end
   end
%Bcs to diagonals

%left
i=1;
for j=1:M
k=(j-1)*N+i;
EE(k)=1;
DD(k)=0;
BB(k)=0;
HH(k)=0;
FF(k)=0;
QQ(k)=phi_L(j);
phi_g(k)=QQ(k);
end
%right
i=N;
for j=1:M
k=(j-1)*N+i;
EE(k)=1;
DD(k)=0;
BB(k)=0;
HH(k)=0;
FF(k)=0;
QQ(k)=phi_R(j);
phi_g(k)=QQ(k);
end
%bottom
j=1;
for i=1:N
k=(j-1)*N+i;
EE(k)=1;
DD(k)=0;
BB(k)=0;
HH(k)=0;
FF(k)=0;
QQ(k)=phi_B(i);
phi_g(k)=QQ(k);
end 
%top
j=M;
for i=1:N
k=(j-1)*N+i;
EE(k)=1;
DD(k)=0;
BB(k)=0;
HH(k)=0;
FF(k)=0;
QQ(k)=phi_T(i);
phi_g(k)=QQ(k);
end


%Looping
count_MSD_41=0;
R2_MSDr_41=1;
R_MSD=zeros(1,neq);

while R2_MSDr_41>Res_Cond && count_MSD_41<20000
%code snippet on page 149 goes here
%step 1
for i=2:N-1
    for j=2:M-1
        k=(j-1)*N+i;
        R_MSD(k)=QQ(k)-EE(k)*phi_g(k)-FF(k)*phi_g(k+1)-HH(k)*phi_g(k+N)-DD(k)*phi_g(k-1)-BB(k)*phi_g(k-N);
    end
end
R2_sum_MSD=0;
for k=1:k
    R2_sum_MSD=R2_sum_MSD+R_MSD(k)*R_MSD(k);
end
R2_MSDr_41=sqrt(R2_sum_MSD);
%step 2
c=zeros(N,M);
for i=2:N-1
    for j=2:M-1
        k=(j-1)*N+i;
        c(k)=EE(k)*R_MSD(k)+FF(k)*R_MSD(k+1)+HH(k)*R_MSD(k+N)+DD(k)*R_MSD(k-1)+BB(k)*R_MSD(k-N);
    end
end
rtc=0;
for k=1:k
rtc=rtc+R_MSD(k)*c(k);
end
alpha=R2_sum_MSD/rtc;

%step 3
phi_g=phi_g+alpha*R_MSD;
count_MSD_41=count_MSD_41+1;
 R2_MSDr_count_41(count_MSD_41)=R2_MSDr_41;
end 

%plotting
%transposing
for i=1:N
    for j=1:M
        k=(j-1)*N+i;
        phi_ij(i,j)=phi_g(k);
    end
end
phi_ij=phi_ij';
for i=1:M/2
t_fact=phi_ij(M+1-i,:);
phi_ij(M+1-i,:)=phi_ij(i,:);
phi_ij(i,:)=t_fact;
end
[X, Y] = meshgrid(x,y(end:-1:1));
%calculate analytical here.
phi_ana=500*exp(-50*((1-X).^2+Y.^2))+100*X.*(1-Y);


figure(1)
contour(X,Y,phi_ij,'ShowText','on','LevelStep',50)
xlabel('x')
ylabel('y')

%R2, error with analyitical contour
figure (2)
semilogy(1:count_MSD_41, R2_MSDr_count_41(1:count_MSD_41))
xlabel('iterations')
ylabel('Residual')

figure (3)
contour(X,Y,phi_ana-phi_ij,'ShowText','on','LevelStep',0.1)
xlabel('x')
ylabel('y')

%% Conjugant Gradient Method

%steps 1-3 before while loop
%step 6 update D vector

N=41;
M=41;
neq=N*M;
phi_g=zeros(1,neq);
Res_cond=10^-6;
%intial variables

del_x=1/(N-1);
del_y=1/(M-1);
Res_Cond=1e-6;
x=linspace(0,1,N);
y=linspace(0,1,M);
%BC's
phi_R=100.*(1-y)+500*exp(-50.*y.^2);
phi_L=500*exp(-50*(1+y.^2));
phi_B=100*x+500*exp(-50*(1-x).^2);
phi_T=500*exp(-50*((1-x).^2+1));
%diagonals from HW4
EE=(2/(del_x^2)+2/(del_y^2))*ones(1,neq);
BB=[zeros(1,N),(-1/del_y^2)*ones(1,neq-N)];
DD=[0,-1/(del_x^2)*ones(1,neq-1)];
FF=[(-1/(del_x^2))*ones(1,neq-1),0];
HH=[(-1/(del_y^2))*ones(1,neq-N),zeros(1,N)];

    for i=1:N
       for j=1:M
       k=(j-1)*N+i;
       QQ(k)=-50000.*exp(-50.*((1-x(i)).^2+y(j).^2)).*(100.*((1-x(i)).^2+y(j).^2)-2);
  
       end
   end
%Bcs to diagonals

%left
i=1;
for j=1:M
k=(j-1)*N+i;
EE(k)=1;
DD(k)=0;
BB(k)=0;
HH(k)=0;
FF(k)=0;
QQ(k)=phi_L(j);
phi_g(k)=QQ(k);
end
%right
i=N;
for j=1:M
k=(j-1)*N+i;
EE(k)=1;
DD(k)=0;
BB(k)=0;
HH(k)=0;
FF(k)=0;
QQ(k)=phi_R(j);
phi_g(k)=QQ(k);
end
%bottom
j=1;
for i=1:N
k=(j-1)*N+i;
EE(k)=1;
DD(k)=0;
BB(k)=0;
HH(k)=0;
FF(k)=0;
QQ(k)=phi_B(i);
phi_g(k)=QQ(k);
end 
%top
j=M;
for i=1:N
k=(j-1)*N+i;
EE(k)=1;
DD(k)=0;
BB(k)=0;
HH(k)=0;
FF(k)=0;
QQ(k)=phi_T(i);
phi_g(k)=QQ(k);
end


%Looping
count_CGM_41=0;
R2_CGMr=1;
R_CGM=zeros(1,neq);
for i=2:N-1
    for j=2:M-1
        k=(j-1)*N+i;
        R_CGM(k)=QQ(k)-EE(k)*phi_g(k)-FF(k)*phi_g(k+1)-HH(k)*phi_g(k+N)-DD(k)*phi_g(k-1)-BB(k)*phi_g(k-N);
    end
end
R2_sum_CGM=0;
for k=1:k
    R2_sum_CGM=R2_sum_CGM+R_CGM(k)*R_CGM(k);
end
%R2_sum_CGM=zeros(1,neq);
D_vect=R_CGM;

while R2_CGMr>Res_Cond && count_CGM_41<20000
% for i=2:N-1
%     for j=2:M-1
%         k=(j-1)*N+i;
%         R_CGM(k)=QQ(k)-EE(k)*phi_g(k)-FF(k)*phi_g(k+1)-HH(k)*phi_g(k+N)-DD(k)*phi_g(k-1)-BB(k)*phi_g(k-N);
%     end
% end

c=zeros(N,M);
for i=2:N-1
    for j=2:M-1
        k=(j-1)*N+i;
        c(k)=EE(k)*D_vect(k)+FF(k)*D_vect(k+1)+HH(k)*D_vect(k+N)+DD(k)*D_vect(k-1)+BB(k)*D_vect(k-N);
    end
end
rtc=0;
for k=1:k
rtc=rtc+D_vect(k)*c(k);
end
alpha=R2_sum_CGM/rtc;

phi_g=phi_g+alpha.*D_vect;

for i=2:N-1
    for j=2:M-1
        k=(j-1)*N+i;
        R_CGM(k)=QQ(k)-EE(k)*phi_g(k)-FF(k)*phi_g(k+1)-HH(k)*phi_g(k+N)-DD(k)*phi_g(k-1)-BB(k)*phi_g(k-N);
    end
end

R2_sum_CGM_2=0;
for k=1:k
    R2_sum_CGM_2=R2_sum_CGM_2+R_CGM(k)*R_CGM(k);
end

beta=R2_sum_CGM_2./R2_sum_CGM;
D_vect=R_CGM+beta*D_vect;


R2_CGMr=sqrt(R2_sum_CGM_2);
R2_sum_CGM=R2_sum_CGM_2;

count_CGM_41=count_CGM_41+1;
 R2_CGMr_count_41(count_CGM_41)=R2_CGMr;
end 


%plotting
%transposing
for i=1:N
    for j=1:M
        k=(j-1)*N+i;
        phi_ij(i,j)=phi_g(k);
    end
end
phi_ij=phi_ij';
for i=1:M/2
t_fact=phi_ij(M+1-i,:);
phi_ij(M+1-i,:)=phi_ij(i,:);
phi_ij(i,:)=t_fact;
end
[X, Y] = meshgrid(x,y(end:-1:1));
%calculate analytical here.
phi_ana=500*exp(-50*((1-X).^2+Y.^2))+100*X.*(1-Y);


figure(4)
contour(X,Y,phi_ij,'ShowText','on','LevelStep',50)
xlabel('x')
ylabel('y')

%R2, error with analyitical contour
figure (5)
semilogy(1:count_CGM_41, R2_CGMr_count_41(1:count_CGM_41))
xlabel('iterations')
ylabel('Residual')

figure (6)
contour(X,Y,phi_ana-phi_ij,'ShowText','on','LevelStep',0.1)
xlabel('x')
ylabel('y')


