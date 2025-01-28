clc
clear all
close all
tic
B=1;
Lambda=[0.8603 3.4256 6.4373 9.5293];

N_x=26;
%minus one becuase 26 nodes but 25 spaces
del_x=1/(N_x-1);
x_vect=0:del_x:1;
%based on stability condition del_t,max=del_x^2/2alpha
del_t_max=del_x^2/2;
del_t=del_t_max/2;
%initialising

phi_p=ones(1,N_x);
phi_n=zeros(1,N_x);
% phi_n(1,1)=del_t/(del_x^2)*phi_n(end)-del_t/(del_x^2)*phi_p(end)+phi_p(end);
% phi_n(1,end)=phi_p(end)*del_t/del_x-(del_t/del_x)*phi_p(end-1)*(1+del_x*B)+phi_p(end);
t_total=0.8;

       a=ones(1,N_x-1)*(-1/(2*del_x^2));
       %a is for j-1 of LHS
       %a(1,1)=-del_t/del_x^2;
       a(1,end)=-1/del_x^2;
       %c is for j+1 terms of LHS
       c=ones(1,N_x)*(-1/(2*del_x^2));
       c(1,1)=-1/del_x^2;
       %d is for phi^n+1 at j of LHS
       d=ones(1,N_x)*(1/del_t+1/del_x^2);
       %d(1,1)=1+1/del_x^2;
       %d(1,end)=1+del_t/del_x*(1+del_x*B);
        d(1,end)=1/del_t+(1+del_x*B)/del_x^2;
RHS=ones(1,N_x);
for i=1:1:t_total/del_t
    for x=2:1:N_x-1
        %iterations
    %     phi_p(x)=phi_p(x)/del_t+0.5*((phi_p(x+1)+phi_p(x-1)-2*phi_p(x))/del_x^2);
    RHS(x)=phi_p(x)/del_t+0.5*((phi_p(x+1)+phi_p(x-1)-2*phi_p(x))/del_x^2);
    end
% RHS(1,1)=del_t/(del_x^2)*phi_n(1)-del_t/(del_x^2)*phi_p(1)+phi_p(1);
% RHS(1,end)=phi_p(end)*del_t/del_x-(del_t/del_x)*phi_p(end-1)*(1+del_x*B)+phi_p(end);
RHS(1,1)=phi_p(1)/del_t-1/del_x^2*phi_p(1)+1/del_x^2*phi_p(2);
RHS(1,end)=phi_p(end)/del_t-((1+del_x*B)/(del_x)^2)*phi_p(end)+phi_p(end-1)/del_x^2;

       

       phi_n=TriDiagS(a,c,N_x,d,RHS);
       
    

         if i*del_t==0.1
             Phi_Num_save_01=phi_n;
         end

         if i*del_t==0.2
             Phi_Num_save_02=phi_n;
         end

         if i*del_t==0.4
             Phi_Num_save_04=phi_n;
         end

         if i*del_t==0.8
                Phi_Num_save_08=phi_n;
         end
         phi_p=phi_n;
       
end

phi_ana=zeros(4,N_x);
t_vect=[0.1 0.2 0.4 0.8];
for i=1:1:4
    %i is for each soln point then next for loop goes through each x
    %location
    for x=1:1:N_x
        for m=1:1:4
            C_n(m)=(4*sin(Lambda(m)))./(2*Lambda(m)+sin(2*Lambda(m)));
            %x_vect not x at end. Because this is the proper it finds the x
            %at the particular time.
            phi_ana(i,x)=phi_ana(i,x)+C_n(m).*exp(-(Lambda(m).^2)*t_vect(i)).*cos(Lambda(m)*x_vect(x));
        end
    end
end
toc
figure (1)
%for t=0.1 0.2 0.4 0.8
plot(x_vect,Phi_Num_save_01)
hold on
grid on
xlabel('x')
ylabel('{\phi}')
plot(x_vect,Phi_Num_save_02)
plot(x_vect,Phi_Num_save_04)
plot(x_vect,Phi_Num_save_08)
hold off
title('Solutions for PDE')
legend('t=0.1','t=0.2','t=0.4','t=0.8')
figure(2)
plot(x_vect,phi_ana(1,:)-Phi_Num_save_01)
xlabel('x')
ylabel('{\phi}_{ana}-{\phi}_{num}')
hold on
grid on
plot(x_vect,phi_ana(2,:)-Phi_Num_save_02)
plot(x_vect,phi_ana(3,:)-Phi_Num_save_04)
plot(x_vect,phi_ana(4,:)-Phi_Num_save_08)
hold off
title('Error Plot')
legend('t=0.1','t=0.2','t=0.4','t=0.8','Location','southwest')
% figure(3)
% plot(x_vect,phi_ana)