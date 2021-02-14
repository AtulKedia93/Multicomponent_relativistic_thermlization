
logEout2 = log10(Eout2);
h = histogram(logEout2,'normalization','pdf', 'BinWidth', KT*0.2);
x=h.BinEdges;y=h.Values; y = log10(10000*y); y = y/sum(y);
x(1)=[];
logRFDM1 = log10(RFDM1); logKE = log10(KE);
figure; hold on;
bar(x,y,'barwidth',0.4)
plot(KE,RFDM1,'r','LineWidth',3);
plot(KE,logRFDM1,'r','LineWidth',3);

axes1 = axes('Parent',figure);
set(axes1,'XTickLabel',{'10^{-3}','10^{-2}','10^{-1}','1','10^{1}','10^{2}','10^{3}'},'YTick',zeros(1,0));

% logRFDM1 = log10(10000*RFDM1); logKE = log10(KE);
% h = histogram(logEout2,'normalization','pdf', 'BinWidth', KT*0.2);
% x=h.BinEdges;y=h.Values; y = log10(10000*y);
% x(1)=[];
% figure; hold on;
% bar(x,y,'barwidth',0.4)
% plot(logKE,logRFDM1,'r','LineWidth',3);