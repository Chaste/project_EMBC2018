% script to calculate steady state ring of chain of springs with Pressure 
% force P radially outwards.

clear 
close all
P=0.5*1.0666e4; %3
k=50; %15
R0 = 0.0015;%1.0;
L = 0.012
n_y=10:1:200;

% Short => y is bigger than L0
% Long => y is smaller than L0 

R_2D_equil_norm = zeros(size(n_y));
R_2D_equil_short = zeros(size(n_y));
R_2D_equil_long = zeros(size(n_y));

n_x=1.3*n_y;
n_z = n_y*12/3/pi/sqrt(3)*2;
for i=1:length(n_y)
    R_2D_equil_norm(i) = fzero(@(R) equil_radius_residual(R,P,k,R0,L,n_z(i),n_y(i)),R0);
    R_2D_equil_short(i) = fzero(@(R) equil_radius_residual(R,P,k,R0,L,n_z(i)/2,n_y(i)),R0);
    R_2D_equil_long(i) = fzero(@(R) equil_radius_residual(R,P,k,R0,L,n_z(i)*2,n_y(i)),R0);
end
 
plot(n_y,1000*R_2D_equil_norm,'k');
hold on
plot(n_y,1000*R_2D_equil_short,'k:','LineWidth',2);
plot(n_y,1000*R_2D_equil_long,'k--');

legend('Equilateral','Stretched','Squashed')


xlabel('n_y')
ylabel('Rquilibrium Radius (mm)')
saveaspngandeps(-1,'EquilibriumRadius',7, 7/5, 9);

% figure 
% plot(L0,y_short,'k',L0,y_norm,'r',L0,y_long,'b'); 
% legend('short','norm','long');
% xlabel('L0');
% ylabel('y');

% Nd = [20,40,80];
% Nz = [15,30,60,];
% 
% R_steady=[0.003020841864943
%          0.001961963139447
%          0.001676058013633
%          0.003
%          0.002
%          0.001
%          0.004
%          0.002
%          0.0015];
% 
% MeshSizes = meshgrid(Nd,ones(3,1))   
% 
% plot(MeshSizes,reshape(R_steady,size(MeshSizes)),'ko')
     

