close all
clc
clear all

tic
figure (999)
% Jacobian
% HW3_81;
% semilogy(1:count_J_81,R2_vect_J_81(1:count_J_81))
hold on
HW3_41;
clearvars -except R2_vectGS_41 R2_vect_J_41 count_J_41 count_GS_41
% Jacobian and Gauss Seidel 41
semilogy(1:count_J_41,R2_vect_J_41(1:count_J_41))
hold on
semilogy(1:count_GS_41,R2_vectGS_41(1:count_GS_41))
% Jacobian and Gauss Seidel 81
HW3_81;
clearvars -except R2_vectGS_81 R2_vect_J_81 count_J_81 count_GS_81 R2_vectGS_41 R2_vect_J_41 count_J_41 count_GS_41
semilogy(1:count_J_81,R2_vect_J_81(1:count_J_81))
semilogy(1:count_GS_81,R2_vectGS_81(1:count_GS_81))
% %ADI and SIP 41
 HW4_41;
 clearvars -except R2_vect_SIP_41 R2_vect_ADI_41 count_ADI_41 count_SIP_41 R2_vectGS_81 R2_vect_J_81 count_J_81 count_GS_81 R2_vectGS_41 R2_vect_J_41 count_J_41 count_GS_41
semilogy(2:2:count_ADI_41,R2_vect_ADI_41(2:2:count_ADI_41))
semilogy(1:count_SIP_41,R2_vect_SIP_41(1:count_SIP_41))
% %ADI and SIP 81
HW4_81;
clearvars -except R2_vect_SIP_81 R2_vect_ADI_81 count_ADI_81 count_SIP_81 R2_vect_SIP_41 R2_vect_ADI_41 count_ADI_41 count_SIP_41 R2_vectGS_81 R2_vect_J_81 count_J_81 count_GS_81 R2_vectGS_41 R2_vect_J_41 count_J_41 count_GS_41
semilogy(1:count_SIP_81,R2_vect_SIP_81(1:count_SIP_81))
semilogy(2:2:count_ADI_81,R2_vect_ADI_81(2:2:count_ADI_81))
hold off
xlabel('iterations')
ylabel('Residual')
legend('Jacobian 41','GS 41','Jacobian 81','GS 81','ADI 41','SIP 41','SIP 81','ADI 81')
toc