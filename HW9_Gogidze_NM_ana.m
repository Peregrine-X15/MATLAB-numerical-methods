clear all
clc
close all
HW3;

del_phi_delx=50000*exp(-50.*((1-x).^2+y.^2)).*(1-x)+100.*(1-y);
del_phi_dely=-50000*exp(-50.*((1-x).^2+y.^2)).*y-100.*x;

QE=del_phi_delx*del_y;
QW=del_phi_delx*del_y;
QS=del_phi_dely*del_x;
QN=del_phi_dely*del_x;

QE=sum(QE);
QW=sum(QW);
QS=sum(QS);
QN=sum(QN);

S_phi= @(x,y) 50000.*exp(-50.*((1-x).^2+y.^2)).*(100.*((1-x).^2+y.^2)-2);

xmin =0; xmax=1;
ymin=0; ymax=1;

P=integral2(S_phi,xmin,xmax,ymin,ymax);

I=QE-QW+QN-QS-P;
