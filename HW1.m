clc
clear all 
close all
%% #1
for a_N=[6 11 21 41]
x=linspace(0,1,a_N);
phi=exp(-x)+(2-exp(-1))*x-1;
figure (1)
plot(x,phi)
hold on
xlabel('x')
ylabel('{\phi}')
grid on
%legend('N=6','N=11','N=21','N=41')
end
hold off
%legend('N=6','N=11','N=21','N=41')

%% 2 use sample code on page 113

%N=[6 11 21 41];
%N=21;
for N=[6 11 21 41]
%N=21;
x_n=linspace(0,1,N);
l=1;
del_x=l/(N-1);
S_phi=exp(-x_n);
S_phi(1,1)=0;
S_phi(1,N)=1;
a_coeff=1./del_x.^2;
d_coeff=-2./del_x.^2;
c_coeff=1./del_x.^2;
% making the tridiagonal matrix 
% coeff_matrix=diag(a_coeff*ones(1,N-1),-1)+diag(d_coeff*ones(1,N))+diag(c_coeff*ones(1,N-1),1);
% coeff_matrix(1,1)=1;
% coeff_matrix(1,2:N)=0;
% coeff_matrix(N,N)=1;
% coeff_matrix(N,N-1)=0;
% vectors of diagonal
a=a_coeff*ones(1,N-1);
a(1,end)=0;
d=d_coeff*ones(1,N);
d(1,1)=1;
d(1,end)=1;
c=c_coeff(ones(1,N-1));
c(1,1)=0;

for i=2:N
xmult=a(i-1)/d(i-1);
d(i)=d(i)-xmult*c(i-1);
S_phi(i)=S_phi(i)-xmult*S_phi(i-1);
end

phi_new=zeros(1,N);
phi_new(1,end)=S_phi(1,end)/d(1,end);

for i=N-1:-1:1
phi_new(i)=(S_phi(i)-c(i)*phi_new(i+1))/d(i);
end
phi=exp(-x_n)+(2-exp(-1))*x_n-1;
% error
E=((phi(2:N-1)-phi_new(2:N-1))./phi(2:N-1))*100;
figure (2)
hold on
plot(x_n(2:N-1),E)
xlabel('x')
ylabel('Error(%)')
grid on

end
hold off
legend('N=6','N=11','N=21','N=41')