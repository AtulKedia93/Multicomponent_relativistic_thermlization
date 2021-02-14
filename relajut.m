syms V v cth E vp;
mcs=.5;
KT=1;
m1=938;
m2=2*m1;
M=m1+m2;

u=(mcs/KT)*(1/besselk(2,mcs/KT));
f=@(E)E.*(E.^2-mcs^2).^.5./(1+exp(E./KT));
n= integral(f,mcs,Inf);
% A=u*exp(-mcs/KT);
v1=@(V,v,cth)(V.^2+2*m2.*V.*v.*cth./M+(m2.*v./M).^2).^.5;
v2=@(V,v,cth)(V.^2-2*m1.*V.*v.*cth./M+(m1.*v./M).^2).^.5;


gp1=@(V,v,cth) 1./sqrt(1-v1(V,v,cth).^2);
gp2=@(V,v,cth) 1./sqrt(1-v2(V,v,cth).^2);
%% max-jutt
% gp1=@(V,v,cth) 1./sqrt(1-v1(V,v,cth).^2);
% g1=@(V,v,cth) 1+0.5*(m1/mcs)*v1(V,v,cth).^2;
% g1=@(V,v,cth) 1+(gp1(V,v,cth)-1).*(m1/mcs); %% relativistic correction gamma-nuclei added
% b1=@(V,v,cth) sqrt(1-1./(g1(V,v,cth).^2));
% f1=@(V,v,cth) (1./(4*pi.*v1(V,v,cth).^2)).* g1(V,v,cth).^2.*b1(V,v,cth).*u.*exp(-g1(V,v,cth).*mcs./KT).*(m1./mcs).*v1(V,v,cth).*gp1(V,v,cth).^3;

% gp2=@(V,v,cth) 1./sqrt(1-v2(V,v,cth).^2);
% g2=@(V,v,cth) 1+0.5*(m2/mcs)*v2(V,v,cth).^2;
% g2=@(V,v,cth) 1+(gp2(V,v,cth)-1).*(m2/mcs); %% relativistic correction gamma-nuclei added
% b2=@(V,v,cth) sqrt(1-1./(g2(V,v,cth).^2));
% f2=@(V,v,cth) (1./(4*pi.*v2(V,v,cth).^2)).* g2(V,v,cth).^2.*b2(V,v,cth).*u.*exp(-g2(V,v,cth).*mcs./KT).*(m2./mcs).*v2(V,v,cth).*gp2(V,v,cth).^3;

%% max boltz 1964
gp=@(vp)1./sqrt(1-vp.^2);
eb1=@(vp) (gp(vp)-1)*m1  + mcs;
pb1=@(vp) sqrt((eb1(vp).^2) - mcs^2);
nb1=@(vp) m1.*vp.*gp(vp).^3.*eb1(vp).*pb1(vp).*(exp(eb1(vp)./KT)+1).^-1;
n11=integral(nb1,0,1);
eb2=@(vp) (gp(vp)-1)*m2  + mcs;

pb2=@(vp) sqrt((eb2(vp).^2) - mcs^2);
nb2=@(vp) m2.*vp.*gp(vp).^3.*eb2(vp).*pb2(vp).*(exp(eb2(vp)./KT)+1).^-1;
n22=integral(nb2,0,1);

% fb=@(vp)nb(vp)/n11;
% 
% 
E1=@(V,v,cth) 0.5*m1.*v1(V,v,cth).^2+mcs;
E2=@(V,v,cth) 0.5*m2.*v2(V,v,cth).^2+mcs;
% E1=@(V,v,cth) (gp1(V,v,cth)-1).*m1+mcs;
% E2=@(V,v,cth) (gp2(V,v,cth)-1).*m2+mcs;
P1=@(V,v,cth) sqrt(E1(V,v,cth).^2-mcs^2);
P2=@(V,v,cth) sqrt(E2(V,v,cth).^2-mcs^2);
% f1=@(V,v,cth) (1/n11)*(1./(4*pi.*v1(V,v,cth).^2)).* m1.*v1(V,v,cth).*E1(V,v,cth).*P1(V,v,cth).*(exp(E1(V,v,cth)./KT)+1).^-1 .*gp1(V,v,cth).^3;
% f2=@(V,v,cth) (1/n22)*(1./(4*pi.*v2(V,v,cth).^2)).* m2.*v2(V,v,cth).*E2(V,v,cth).*P2(V,v,cth).*(exp(E2(V,v,cth)./KT)+1).^-1 .*gp2(V,v,cth).^3;

f1=@(V,v,cth) (1/n11)*(1./(4*pi.*v1(V,v,cth).^2)).* m1.*v1(V,v,cth).*E1(V,v,cth).*P1(V,v,cth).*(exp(E1(V,v,cth)./KT)+1).^-1 ;
f2=@(V,v,cth) (1/n22)*(1./(4*pi.*v2(V,v,cth).^2)).* m2.*v2(V,v,cth).*E2(V,v,cth).*P2(V,v,cth).*(exp(E2(V,v,cth)./KT)+1).^-1 ;

%% final rel-vel calc
f=@(v)vpaintegral(vpaintegral(2*pi.*V.^2*f2(V,v,cth).*f1(V,v,cth),cth,[-1 1]),V,[0 1]);
% fun=@(V,v,cth)2*pi.*V.^2.*f2(V,v,cth).*f1(V,v,cth);
% f=integral3(fun,0,1,0,1,-1,1);

parfor i=1:200 %% problem for relativistic gamma corrected nuclei for high values of "i"
   x(i)=(i-1)*0.005;
   y(i)=double(4*pi*(x(i))^2*f(x(i)));
end

mu=(m1*m2/(m1+m2));
mb=@(l)4*pi.*(l).^2 .* (mu/(2*pi*KT))^1.5.*exp(-mu.*(l).^2/(2*KT));

width=7;
height=5;
alw=0.75;
fsz=12;
lw=1.5;
msz=6;
figure(1);
pos=get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) width*100 height*100]);
% mb=@(l)4*pi.*(l).^2 .* (1/n11).*exp(-mu.*(l).^2/(2*KT));
plot(x,mb(x),'r','linewidth',lw,'markersize',msz);
hold on
plot(x,y,'b','linewidth',lw,'markersize',msz);
set(gca,'fontsize',fsz,'linewidth',alw);
ax=legend('Standard Maxwellian','Distorted');
%leg=findobj(ax,'type','text');
ax.FontSize=fsz;
legend boxoff;
ylabel('f(v)','fontsize',fsz)
xlabel('velocity','fontsize',fsz)
text(0.18,12.5,'KT = 1.0 MeV','fontsize',fsz)
xlim([0 0.3])
set(gcf,'inverthardcopy','on');
set(gcf,'paperunits','inches');
papersize=get(gcf,'papersize');
left=(papersize(1)-width)/2;
bottom=(papersize(2)-height)/2;
myfiguresize=[left,bottom,width,height];
set(gcf,'paperposition',myfiguresize);
print('1p0relvel','-dpng','-r300')
savefig('1p0relvel')