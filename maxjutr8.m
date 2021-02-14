i=1;
zi=1;
zj=1;
m1=1.007276466879;      % incoming masses, proton
m2=2.01410178;          % incoming masses, deuterium
mu=m1*m2/(m1+m2);        % reduced mass
K=8.6173303e-2;
b=0.9898*zi*zj*mu^0.5;   % factors of reaction rate/ equation (3)
s0=-0.1857;                 % E to be entered in MeV
s1=8.675;
s2=1.922;
s3=-2.197;
c=(3*10^10)*10^-24;
mcs=0.5;
mcs2=931.5*mu;
q=zeros(10000,1);

fileID = fopen('f20.txt','w');

for T=.001:.001:10
    u=exp(-mcs/(K*T))*(c/(K*T))*(1/mcs)^2*(1/besselk(2,(mcs/(K*T))));
    g=@(E) sqrt(2*E./mcs2).*(1./E).*(1-(mcs./(E+mcs)).^2).^0.5;
    ex=@(E) exp(-(E./(K*T))-(b./sqrt(E)));
    S=@(E) (s0+s1*E+s2*E.^2+s3*E.^3);
    fun=@(E)u.*g(E).*ex(E).*S(E);
    q(i) = integral(fun,mcs,Inf);
    if(q(i)<1.0e-25||isnan(q(i)))
        q(i)=0.0;
    end
    fprintf(fileID,'%f %f\n',T, q(i));
    i=i+1;
end
fclose(fileID);
j=linspace(.001,10,10000);
loglog(j,q);