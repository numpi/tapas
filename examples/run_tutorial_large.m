minm = 5;
maxm = 15;
results = zeros(maxm-minm,3);
for m = minm : maxm
    [MTTF, CMTTF, CMTTB] = tutorial_large(m);
    results(m-minm+1,1) = MTTF;
    results(m-minm+1,2) = CMTTF;
    results(m-minm+1,3) = CMTTB;
end

slg=semilogy(minm:maxm,results(:,1),'--',minm:maxm,results(:,2),minm:maxm,results(:,3),'o');
xlabel('m (number of copies of A1 in Figure 1)');
ylabel('hours');
legend('MTTF','CMTTF','CMTTB','location','northwest');
slg(1).Color = [0 0 0];
slg(2).Color = [0 0 0];
slg(3).Color = [0 0 0];