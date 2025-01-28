clear all
clc
close all
tic
N_x=26;
%follow lecture 19, make additional phi. Initial condition is first phi
%then updated. Consult slides 7 and 10 
%issue might be with storing
del_x=1./(N_x-1);
% part a
del_t_max=del_x;
del_t=del_t_max*2;
x_vect=0:del_x:1;
t_total=1.2;
phi_p=zeros(1,N_x);
%initial condition
phi_c=ones(1,N_x).*(sin(pi*x_vect));
%*(sin(pi*x_vect)+del_t*0.25*sin(2*pi*x_vect)+0.5*del_t^2/del_x^2);
phi_n=zeros(1,N_x);
for x=2:1:N_x-1
        phi_n(x)=phi_c(x)+del_t*0.25*sin(2*pi*x_vect(x))+0.5*del_t^2/del_x^2*(phi_c(x+1)-2*phi_c(x)+phi_c(x-1));
end
           phi_p=phi_c;
          phi_c=phi_n;
for i=2:0.1:t_total/del_t%, doesnt equal 3000 for some reason but it does in console
for j=2:1:N_x-1
    phi_n(j)=2*phi_c(j)-phi_p(j)+(del_t^2/del_x^2)*(phi_c(j-1)-2*phi_c(j)+phi_c(j+1));
end
    phi_n(1)=0;
    phi_n(end)=0;
         if i*del_t==0.4
             Phi_Num_save_04_b=phi_n;
         end

         if i*del_t==0.8
             Phi_Num_save_08_b=phi_n;
         end

         if i*del_t==1.2
             Phi_Num_save_12_b=phi_n;
         end
          phi_p=phi_c;
          phi_c=phi_n;
   
end
toc
figure (1)
plot(x_vect,Phi_Num_save_04_b)
xlabel('x')
ylabel('{\phi}')
hold on
grid on
plot(x_vect,Phi_Num_save_08_b)
plot(x_vect,Phi_Num_save_12_b)
hold off
title('Solutions for PDE')
legend('t=0.4','t=0.8','t=1.2','Location','southeast')