% clear all
% clc
% close all
tic
N_x=26;
%follow lecture 19, make additional phi. Initial condition is first phi
%then updated. Consult slides 7 and 10 
del_x=1/(N_x-1);
% part a
del_t_max=del_x;
del_t=del_t_max/2;
x_vect=0:del_x:1;
t_total=1.2;
phi_p=zeros(1,N_x);
%initial condition
phi_c=ones(1,N_x).*(sin(pi*x_vect));
phi_n=zeros(1,N_x);
%setup TDMA here
a=ones(1,N_x-1).*(-del_t^2./del_x^2);
d=ones(1,N_x).*(1+del_t^2./del_x^2);
c=ones(1,N_x-1).*(-del_t^2./del_x^2);
c(1,1)=0;
a(1,end)=0;
d(1,end)=1;
d(1,1)=1;
b=sin(pi*x_vect)+del_t*0.25*sin(2*pi*x_vect);
b(1,1)=0;
b(1,end)=0;
%for x=2:1:N_x-1
phi_n=TriDiagS(a,c,N_x,d,b);
%end
phi_p=phi_c;
phi_c=phi_n;
%TDMA for the next equation using results from above.
for i=2.0:1.0:3000 %t_total/del_t, doesnt equal 3000 for some reason but it does in console
%TDMA for second eqn goes here.
    a_2=ones(1,N_x-1)*(-del_t^2/del_x^2);
    a_2(1,end)=0;
    d_2=ones(1,N_x)*(1+2*del_t^2/del_x^2);
    c_2=ones(1,N_x-1)*(-del_t^2/del_x^2);
    d_2(1,1)=1;
    d_2(1,end)=1;
    c_2(1,1)=0;
    b_2=2*phi_c-phi_p;
    b_2(1,1)=0;
    b_2(1,end)=0;

    phi_n=TriDiagS(a_2,c_2,N_x,d_2,b_2);
         if i*del_t==0.4
             Phi_Num_save_04=phi_n;
         end

         if i*del_t==0.8
             Phi_Num_save_08=phi_n;
         end

         if i*del_t==1.2
             Phi_Num_save_12=phi_n;
         end
          phi_p=phi_c;
          phi_c=phi_n;
end
toc
% figure (1)
% plot(x_vect,Phi_Num_save_04)
% xlabel('x')
% ylabel('{\phi}')
% hold on
% grid on
% plot(x_vect,Phi_Num_save_08)
% plot(x_vect,Phi_Num_save_12)
% hold off
% title('Solutions for PDE')
% legend('t=0.4','t=0.8','t=1.2','Location','southeast')