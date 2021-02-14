syms V v cth E;
mcs=.5;
KT=1;
m1=938;
m2=2*m1;
M=m1+m2;

u=(1/KT)*(1/besselk(2,mcs/KT));
A=u*exp(-mcs/KT);
k1=@(V,v,cth) 0.5*m1*(V.^2+2*m2.*V.*v.*cth./M+(m2.*v./M).^2);
k2=@(V,v,cth) 0.5*m2*(V.^2-2*m1.*V.*v.*cth./M+(m1.*v./M).^2);
x1=@(V,v,cth) ((k1(V,v,cth)+mcs)./mcs).^2;
x2=@(V,v,cth) ((k2(V,v,cth)+mcs)./mcs).^2;
v1=@(V,v,cth) sqrt(k1(V,v,cth)./(m1*.5));
v2=@(V,v,cth) sqrt(k2(V,v,cth)./(m2*.5));

f=@(E)E.*(E.^2-mcs^2)^.5./(1+exp(E./KT));
n1= vpaintegral(f(E),E,[mcs inf]);
%maxwell-bolzmann from 1964
f1=@(V,v,cth) (m1./(4*pi.*v1(V,v,cth))).*(1/n1).*(.5*m1.*v1(V,v,cth).^2+mcs).*((.5*m1.*v1(V,v,cth).^2+mcs).^2-mcs^2)^.5./(1+exp((.5*m1.*v1(V,v,cth).^2+mcs)./KT));
f2=@(V,v,cth) (m2./(4*pi.*v2(V,v,cth))).*(1/n1).*(.5*m2.*v2(V,v,cth).^2+mcs).*((.5*m2.*v2(V,v,cth).^2+mcs).^2-mcs^2)^.5./(1+exp((.5*m2.*v2(V,v,cth).^2+mcs)./KT));

%maxwell Juttner
% f1=@(V,v,cth)   ( m1./(4*pi.*v1(V,v,cth)) )   .*A.* sqrt((x1(V,v,cth)-1).*x1(V,v,cth))  .* exp(-k1(V,v,cth)./KT);
% f2=@(V,v,cth)   ( m2./(4*pi.*v2(V,v,cth)) )   .*A.* sqrt((x2(V,v,cth)-1).*x2(V,v,cth))  .* exp(-k2(V,v,cth)./KT);

% k1=.5*m1*(V^2+2*m2*V*v*cth/M+(m2*v/M)^2);
% k2=.5*m2*(V^2-2*m1*V*v*cth/M+(m1*v/M)^2);
% f1=@(V,v,cth)A*(m1./(4*pi*sqrt(V^2+2*m2*V*v*cth/M+(m2*v/M)^2))).*(V.^2+2*m2.*V.*v.*cth/M+(m2.*v./M).^2).^.5.*(.5.*m1.*(V.^2+2*m2.*V.*v.*cth./M+(m2.*v/M).^2)./mcs+1).*((.5.*m1.*(V.^2+2.*m2.*V.*v.*cth./M+(m2.*v./M).^2)./mcs).^2+(2.*.5*m1.*(V.^2+2*m2.*V.*v.*cth./M+(m2.*v./M).^2)./mcs)).^.5.*exp(-.5*m2*(V.^2+2*m2.*V.*v.*cth./M+(m2.*v./M).^2)./KT);
% f2=@(V,v,cth)A*(m2./(4*pi*sqrt(V^2-2*m1*V*v*cth/M+(m1*v/M)^2))).*(V.^2-2*m1.*V.*v.*cth/M+(m1.*v./M).^2).^.5.*(.5.*m2.*(V.^2-2*m1.*V.*v.*cth./M+(m1.*v/M).^2)./mcs+1).*((.5.*m2.*(V.^2-2.*m1.*V.*v.*cth./M+(m1.*v./M).^2)./mcs).^2+(2.*.5*m2.*(V.^2-2*m1.*V.*v.*cth./M+(m1.*v./M).^2)./mcs)).^.5.*exp(-.5*m2*(V.^2-2*m1.*V.*v.*cth./M+(m1.*v./M).^2)./KT);
% f11=@(V,v)int(f1(V,v,cth).*f2(V,v,cth),cth,-1,1);
f=@(v)vpaintegral(vpaintegral(2*pi.*V.^2*f2(V,v,cth).*f1(V,v,cth),cth,[-1 1]),V,[0 1]);
% fin=vpaintegral(4*pi*v.^2.*f(v),v,[0 1]);
parfor i=1:101
   x(i)=(i-1)*0.005;
   y(i)=double(4*pi*((i-1)*0.005)^2*f((i-1)*0.005));
end
mu=(m1*m2/(m1+m2));
mb=@(l)4*pi.*(l).^2 .* (mu/(2*pi*KT))^1.5.*exp(-mu.*(l).^2/(2*KT));
plot(x,mb(x),'r');
% hold on
% mj=@(z) mu.*z.^2./(1-z.^2).^2.5 * (1/KT)* (1/besselk(2,mu/KT)).*exp(-(mu./sqrt(1-z.^2))./KT);
% plot(x,mj(x),'b');
hold on
plot(x,y,'g');