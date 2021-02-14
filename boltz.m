syms p
mcs=.51;
KT=0.1;
m1=939;
i=1;
% for KT=.01:.01:.1
fd=@(E)(E+mcs).*(E.*(E+2*mcs)).^0.5./(exp(mcs/KT)*exp(E./KT)+1);
n= integral(fd,0,Inf);
fd1=@(E)(1/n)*(E+mcs).*(E.*(E+2*mcs)).^0.5./(exp(mcs/KT)*exp(E./KT)+1);
avk=@(E)E.*(1/n).*(E+mcs).*(E.*(E+2*mcs)).^0.5./(exp(mcs/KT)*exp(E./KT)+1);

% fdd=@(E)(E+mcs).*(E.*(E+2*mcs)).^0.5./(exp(mcs/KT)*exp(E./KT));
% nd= integral(fdd,0,Inf);
% fdd1=@(E)(1/nd)*(E+mcs).*(E.*(E+2*mcs)).^0.5./(exp(mcs/KT)*exp(E./KT));
% avk1=@(E)E.*(1/nd).*(E+mcs).*(E.*(E+2*mcs)).^0.5./(exp(mcs/KT)*exp(E./KT));

av=integral(avk,0,inf);
% some(i)=(3/2)*KT;
% i=i+1;
% end
% z=mcs/KT;
pre=av/KT;
% % pre=z*besselk(3,z)/besselk(2,z)-1-z;
% 
gam=(2/3)*pre;

Fb=@(E)sqrt(E).*(1./(gam*KT)).^1.5.*exp(-E./(gam.*KT));
nn=integral(Fb,0,inf);
Fb1=@(E)sqrt(E).*(1./(gam*KT)).^1.5.*exp(-E./(gam.*KT))*(1/nn);
%Fbb=@(E)sqrt(E).*(1./(gam*KT)).^1.5./(exp(E./(gam.*KT))+1);
%nnn=integral(Fbb,0,inf);
%Fbb1=@(E)sqrt(E).*(1./(gam*KT)).^1.5./(exp(E./(gam.*KT))+1)*(1/nnn);
fplot(fd1,[0 1])
hold on
%fplot(Fbb1,[0 1])
fplot(Fb1,[0 1],'--')