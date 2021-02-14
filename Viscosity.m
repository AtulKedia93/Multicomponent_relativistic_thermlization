% theta = -pi:0.01:pi;
% a = 0.05;
% b = 0.7;
% phi = atan2(b*sin(theta),(a+b*cos(theta)));
% f_phi = (a^2+b^2+2*a*b*cos(theta))./(b^2+a*b*cos(theta));
% % figure;plot(theta,phi)
% % figure;plot(phi,f_phi)
% y=theta;
% f2=-log(-(a - a*exp(abs(y)*2i) + (a^2 + a^2*exp(abs(y)*4i) + 4*b^2*exp(abs(y)*2i) - 2*a^2*exp(abs(y)*2i))^(1/2))/(2*b))*(dirac(abs(y)) - dirac(abs(y) - pi/2))*1i + log((- a + a*exp(abs(y)*2i) + (a^2 + a^2*exp(abs(y)*4i) + 4*b^2*exp(abs(y)*2i) - 2*a^2*exp(abs(y)*2i))^(1/2))/(2*b))*(dirac(abs(y) - pi) - dirac(abs(y) - pi/2))*1i + ((heaviside(abs(y)) - heaviside(abs(y) - pi/2))*(a*exp(abs(y)*2i)*2i - (- a^2*exp(abs(y)*2i)*4i + a^2*exp(abs(y)*4i)*4i + b^2*exp(abs(y)*2i)*8i)/(2*(a^2 + a^2*exp(abs(y)*4i) + 4*b^2*exp(abs(y)*2i) - 2*a^2*exp(abs(y)*2i))^(1/2)))*1i)/(a - a*exp(abs(y)*2i) + (a^2 + a^2*exp(abs(y)*4i) + 4*b^2*exp(abs(y)*2i) - 2*a^2*exp(abs(y)*2i))^(1/2)) + ((heaviside(abs(y) - pi) - heaviside(abs(y) - pi/2))*(a*exp(abs(y)*2i)*2i + (- a^2*exp(abs(y)*2i)*4i + a^2*exp(abs(y)*4i)*4i + b^2*exp(abs(y)*2i)*8i)/(2*(a^2 + a^2*exp(abs(y)*4i) + 4*b^2*exp(abs(y)*2i) - 2*a^2*exp(abs(y)*2i))^(1/2)))*1i)/(- a + a*exp(abs(y)*2i) + (a^2 + a^2*exp(abs(y)*4i) + 4*b^2*exp(abs(y)*2i) - 2*a^2*exp(abs(y)*2i))^(1/2))
% plot(y,f2)

% your solution f(phi) in terms of (theta)
% clear
% a=0.05;
% b=0.0501;
% theta = -pi:0.0001:pi;
% phi = atan2(b*sin(theta),(b*cos(theta)-a));
% %theta = mod(atan2(b*sin(theta),(b*cos(theta)-a)),2*pi);
% f_phi = (a^2+b^2-2*a*b*cos(theta))./(b^2-a*b*cos(theta));
% theta1=randpdf(f_phi,phi,[1,1]);
% plot(phi,f_phi)
% hold;

% Nishanth's solutions 5/2/19
% theta = 0:0.0001:2*pi;
% for i=1:length(theta)
%     if(a+b*cos(theta(i))>0)
%         phi(i) = atan(b*sin(theta(i))./(a+b*cos(theta(i))));
%         if(phi(i)<0)
%             phi(i)=2*pi+phi(i);
%         end
%     else
%         phi(i) = pi-atan(b*sin(pi-theta(i))./(-a+b*cos(pi-theta(i))));
%         phi(i) = pi+atan(b*sin(theta(i))./(a+b*cos(theta(i))));
%         if(phi(i)<0)
%             phi(i)=2*pi+phi(i);
%         end
%     end
% end
% c=a/b;
% f_phi=(c^2+2*c*cos(theta)+1)./(1+c*cos(theta));
% plot(phi,f_phi,'*')
% hold

% syms x y a b
% a=0.04;
% b=0.06;
% c=a/b;
% v=-pi:0.01:pi;
% eqn = tan(y)==b*sin(x)/(b*cos(x)+a);
% g=solve(eqn,x);
% if length(g)>1
%     g1(y)=g(1);
%     g2(y)=g(2);
% else
%     g1(y)=g;
%     g2(y)=g;
% end
% gphi(y)=(g2(y))*(heaviside(y)-heaviside(y-pi/2))+(g1(y))*(heaviside(y-pi/2)-heaviside(y-pi));
% gne(y)=diff(gphi,y);
% gne2(y)=gne(abs(y));
% x=gne2(v);
% v=v+pi;

% a=0.04;
% b=0.06;
% phi=-3*pi/4;
% count_phi=0;
% syms theta1
% eqn = tan(phi) == b*sin(theta1)/(b*cos(theta1)-a);
% sol = double(abs(solve(eqn,theta1)));
% if tan(phi) == b*sin(sol(1))/(b*cos(sol(1))-a)
%     theta1 = sol(1);
% elseif tan(phi) == b*sin(sol(2))/(b*cos(sol(2))-a)
%     theta1 = sol(2);
% else
%     count_phi=count_phi+1;
% end

% a = 0.03;
% b = 0.06;
% G = 1/sqrt(1-a^2);
% theta = 1;
% phi = atan2(b*sin(theta),G*(b*cos(theta)-a));
% theta1 = atan2(b*sin(phi),G*(b*cos(phi)+a));

% syms a b phi theta
% eqn =0;
% sol = 0;
% sol2 = 1;
% parfor i=1:4
% eqn = tan(phi) == b*sin(theta)/(b*cos(theta)-a);
% sol = abs(solve(eqn,theta))
% sol2=double(subs(subs(subs(sol,b,0.5),a,0.3),phi,0.5))
% end
% delete(gcp)

% parfor ii = 1:4
% x = rand(10,10);
% y = ones(1,3);
% parsave(sprintf('output%d.mat', ii), x, y);
% end
% delete(gcp)