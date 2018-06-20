function res = equil_radius_residual(R,P,k,R0,L,n_z,n_y)

sinpin = sin(pi/n_y);
cospin = cos(pi/n_y);
y = 2*R*sinpin;
y0 = 2*R0*sinpin;
%x = sqrt(3)/2* y0;
x = L/(n_z -1);


xhat = x;
yhat = R*sinpin;
zhat = R*(1-cospin);                    

yhat0 = R0*sinpin;
zhat0 = R0*(1-cospin);

l_hat = sqrt(xhat*xhat + yhat*yhat + zhat*zhat);
l0_hat = sqrt(xhat*xhat + yhat0*yhat0 + zhat0*zhat0);

res = P*x*y ...
     - 4*k*sinpin*sinpin*(R-R0) ...
     - 4*k*zhat*(1-l0_hat/l_hat);
