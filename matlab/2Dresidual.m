function re = 2Dresidual(R)

R0= 0.0015;
P=1.0666e4;
k=100;
n=60;

res = 8*R^2*(1-sqrt(2)*R0/(sqrt(R0^2+R^2))*(1-(cos(pi/n))^2)+(1-P*R0/k)*R - R0;