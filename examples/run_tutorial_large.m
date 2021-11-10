minm = 5;
maxm = 15;
% results = zeros(maxm-minm,3);
% for m = minm : maxm
%     [MTTF, CMTTF, CMTTB] = tutorial_large(m);
%     results(m-minm+1,1) = MTTF;
%     results(m-minm+1,2) = CMTTF;
%     results(m-minm+1,3) = CMTTB;
% end

slg=semilogy(minm+2:maxm+2,results(:,1),'--',minm+2:maxm+2,results(:,2),minm+2:maxm+2,results(:,3),'o');
xticks(minm+2:maxm+2);
axis([minm+2,maxm+2,1.0e+1,1.5e+5])
xlabel('n (number of submodels)');
ylabel('hours');
legend('MTTF','CMTTF','CMTTB','location','northwest');
slg(1).Color = [0 0 0];
slg(2).Color = [0 0 0];
slg(3).Color = [0 0 0];