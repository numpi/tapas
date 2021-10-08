mink = 5;
maxk = 15;
results = zeros(maxk-mink,3);
for k = mink : maxk
    [MTTF, CMTTF, CMTTB] = tutorial_large(k);
    results(k-mink+1,1) = MTTF;
    results(k-mink+1,2) = CMTTF;
    results(k-mink+1,3) = CMTTB;
end

xlabel('k');
ylabel('hours');
plot(5:15,results(:,1));