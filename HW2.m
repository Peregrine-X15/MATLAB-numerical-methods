clear all
clc
close all
tic
%arbitrary c value in exp
cc=2;

N=41;
%phi_g=linspace(0,1,N);
%phi_g(1,1)=0;
%phi_g(1,N)=1;
phi_g=linspace(0,1,N);
% phi=zeros(1,N);
% phi(1,1)=0;
% phi(1,N)=1;
x_n=linspace(0,1,N);


l=1;
del_x=l/(N-1);
a_coeff=1./del_x^2;
d_coeff=-2./del_x^2;
c_coeff=1./del_x^2;



a=a_coeff*ones(1,N-1);
a(1,end)=0;
d=d_coeff*ones(1,N);
d(1,1)=1;
d(1,end)=1;
c=c_coeff*(ones(1,N-1));
c(1,1)=0;
count=0;
R2=1;
while R2>10^-10 && count<1000
S_phi=exp(cc*phi_g);
S_phi(1,1)=0;
S_phi(1,N)=1;
phi_star=TriDiagS(a,c,N,d,S_phi);
phi_g=phi_star;
for n=2:40
RHS=exp(cc*phi_star(n));
LHS=(phi_star(n+1)-2*phi_star(n)+phi_star(n-1))/del_x^2;
Res_diff=RHS-LHS;
end
R2=sqrt(sum(Res_diff.^2));
count=count+1;
R2_vect(count)=R2;
end
figure (1)
semilogy(1:count,R2_vect)
xlabel('iterations')
ylabel('Residual')

figure (2)
plot(x_n,phi_star)
xlabel('x')
ylabel('{\phi}')
toc
% filename ='unLinearized data.xlsx';
% ansa='y';
% if ansa=='y'
% xlswrite([filename],[(1:count)' R2_vect(1:count)'],'A');
% xlswrite([filename],[x_n' phi_star'],'C');
% end
%% Linearised source term
tic
clear all
clc
cc=2;
N=41;
phi_g=linspace(0,1,N);

x_n=linspace(0,1,N);


l=1;
del_x=l/(N-1);
a_coeff=-1./del_x^2;
d_coeff=2./del_x^2;
c_coeff=-1./del_x^2;



a=a_coeff*ones(1,N-1);
a(1,end)=0;
d=d_coeff*ones(1,N);
d(1,1)=1;
d(1,end)=1;
c=c_coeff*(ones(1,N-1));
c(1,1)=0;
count=0;
R2=1;

while R2>10^-10 && count<1000

d=d_coeff*ones(1,N);
d(1,1)=1;
d(1,end)=1;
%step 4 of linearisation
Sp=-cc*exp(cc.*phi_g);
Sc=-exp(cc.*phi_g)+cc.*phi_g.*exp(cc.*phi_g);
for ii=2:N-1
 %Sp=-cc*exp(cc.*phi_g(ii));
% Sc=-exp(-cc.*phi_g(ii))+cc.*phi_g(ii).*exp(cc.*phi_g(ii));
d(ii)=d(ii)-min([0,Sp(ii)]); 
 Sc(ii)=Sc(ii)+max([0,Sp(ii)])*phi_g(ii);
end
S_phi=Sc;
S_phi(1,1)=0;
S_phi(1,N)=1;
phi_star=TriDiagS(a,c,N,d,S_phi);
phi_g=phi_star;

for n=2:N-1
%RHS=-exp(cc*phi_star(n))*phi_star(n-1)-exp(cc*phi_star(n))+cc*phi_star(n)*exp(cc*phi_star(n));
RHS=-exp(cc*phi_star(n));
LHS=(2*phi_star(n)-phi_star(n+1)-phi_star(n-1))/del_x^2;
Res_diff=RHS-LHS;
end

R2=sqrt(sum(Res_diff.^2));
count=count+1;
R2_vect(count)=R2;
end

figure (3)
semilogy(1:count,R2_vect(1:count))
xlabel('iterations')
ylabel('Residual')

figure (4) 
plot(x_n,phi_star)
xlabel('x')
ylabel('{\phi}')
toc
% filename ='Linearized data.xlsx';
% ansa='y';
% if ansa=='y'
% xlswrite([filename],[(1:count)' R2_vect(1:count)'],'A');
% xlswrite([filename],[x_n' phi_star'],'C');
% end
% filename ='Linearized data.xlsx';
% writematrix(R2_vect(1:count),filename,'Sheet',1)
% writematrix((1:count),filename,'Sheet',1)
% writematrix(x_n,filename,'Sheet',1)
% writematrix(phi_star,filename,'Sheet',1)
%% Newton's Method
tic
clear all 
clc

cc=2;
N=41;
phi_g=linspace(0,1,N);

x_n=linspace(0,1,N);
l=1;
del_x=l/(N-1);

deri_im=(1/del_x^2)*ones(1,N-1);
deri_im(1,end)=0;
deri_ip=(1/del_x^2)*ones(1,N-1);
deri_ip(1,1)=0;

count=0;
R2=1;
f_phi=zeros(1,N);
while R2>10^-10 && count<1000

for p=2:N-1
f_phi(p)=(phi_g(p+1)-2*phi_g(p)+phi_g(p-1))./del_x^2-exp(cc.*phi_g(p));

deri_i(p)=(-2/del_x^2)-cc*exp(cc*phi_g(p));
end
deri_i(1,1)=1;
deri_i(1,N)=1;

del_phi=TriDiagS(deri_im,deri_ip,N,deri_i,-f_phi);
phi_star=phi_g+del_phi;
phi_g=phi_star;
R2=sqrt(sum(f_phi.^2));
count=count+1;
R2_vect(count)=R2;
end

figure (5)
semilogy(1:count,R2_vect(1:count))
xlabel('iterations')
ylabel('Residual')

figure (6)
plot(x_n,phi_star)
xlabel('x')
ylabel('{\phi}')

% ansa='y';
% filename ='Newtons_Method.xlsx';
% if ansa=='y'
% xlswrite([filename],[(1:count)' R2_vect(1:count)'],'D');
% xlswrite([filename],[x_n' phi_star'],'E');
% end
toc