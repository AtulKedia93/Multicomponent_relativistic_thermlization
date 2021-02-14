syms E;
KT=.01;
mcs=.5;
n=@(E) (E+mcs).*sqrt((E+mcs).^2-mcs^2).*(exp((E+mcs)./KT)+1).^-1;
n1=integral(n,0,Inf);
f=@(E)(1/n1)*n(E);
x=linspace(0,10*KT,100);

% mj=@(E) (1/KT)*(1/besselk(2,mcs/KT))*(E./mcs).^2 .*sqrt(1-1./(E./mcs).^2).*exp(-E./KT);
mb=@(E)2.*sqrt(E./pi).*(1/KT)^1.5.*exp(-E./KT);

plot(x,f(x),'--g','MarkerSize',10,'Linewidth',1)
hold on
plot(x,mb(x),'.r','MarkerSize',10,'Linewidth',.5)
legend('Relativistic','Standard Maxwellian')
ylabel('f(KE)')
xlabel('Kinetic Energy (MeV)')