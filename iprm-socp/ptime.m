clc; clear;
fileID1 = fopen('./time.txt','r');

A_t = fscanf(fileID1,'%f');

numsolver = 2; 
num_of_probs = size(A_t,1)/numsolver;
disp(size(A_t,1))
disp(num_of_probs)
c_t = zeros(num_of_probs,numsolver);

for i = 1:numsolver
    c_t(:, i) = A_t(1 + num_of_probs * (i - 1):num_of_probs * i, 1);
% c_t(:,1) = A_t(1:num_of_probs,1);      %2 norm
% c_t(:,2) = A_t(num_of_probs+1:2*num_of_probs,1);   
% c_t(:,3) = A_t(2*num_of_probs+1:3*num_of_probs,1); 
% c_t(:,4) = A_t(3*num_of_probs+1:4*num_of_probs,1);   
%c_t(:,5) = A_t(4*num_of_probs+1:end,1); 
end


c_min_t = zeros(num_of_probs,1);

for i = 1:num_of_probs
    c_min_t(i) = min(c_t(i,:));
end

r_t = zeros(num_of_probs,numsolver);

for j = 1:numsolver
    for i = 1:num_of_probs
        r_t(i,j) = c_t(i,j)/c_min_t(i);
       % r_it(i,j) = c_it(i,j)/c_min_it(i);
       % r_itall(i,j) = c_itall(i,j)/c_min_itall(i);
    end
end
N = 10;
t_t = linspace(1,N,100);
perf_func_t = zeros(100,numsolver);
for k = 1:100
    for j = 1:numsolver
        perf_func_t(k,j) = sum(r_t(:,j) <= t_t(k))/num_of_probs;
    end
end
%plot(t_t,perf_func_t(:,1),'m--',t_t,perf_func_t(:,2),'b.-',t_t,perf_func_t(:,3),'k:','LineWidth',1.5);
%plot(t_t,perf_func_t(:,1),'m--',t_t,perf_func_t(:,2),'b.-','LineWidth',1.5);
%plot(t_t,perf_func_t(:,1),'m--',t_t,perf_func_t(:,2),'b.-',t_t,perf_func_t(:,3),'k:',t_t,perf_func_t(:,4),'c-',t_t,perf_func_t(:,5),'r-.','LineWidth',1.5);
%semilogx(t_t,perf_func_t(:,1),'m--',t_t,perf_func_t(:,2),'b.-',t_t,perf_func_t(:,3),'k:',t_t,perf_func_t(:,5),'c-',t_t,perf_func_t(:,4),'r-.','LineWidth',1.5);
%semilogx(t_t,perf_func_t(:,1),'m--',t_t,perf_func_t(:,2),'b.-',t_t,perf_func_t(:,3),'k:',t_t,perf_func_t(:,4),'c-','LineWidth',1.5);
semilogx(t_t,perf_func_t(:,1),'m--',t_t,perf_func_t(:,2),'b.-','LineWidth',1.5);
%plot(t_t,perf_func_t(:,1),'m--',t_t,perf_func_t(:,2),'b.-',t_t,perf_func_t(:,3),'k:',t_t,perf_func_t(:,4),'c-','LineWidth',1.5);
%semilogx(t_t,perf_func_t(:,1),'m--',t_t,perf_func_t(:,2),'k:','LineWidth',1.5);
 % plot(t_t,perf_func_t(:,1),'m--',t_t,perf_func_t(:,2),'b.-');
axis([1 N 0 1 ]);
%lgd = legend('IP-PMM','IPRM','OSQP','Location','southeast');
%lgd = legend('IP-PMM','IPRM','Location','southeast');
%lgd = legend('iIPRM','IP-PMM','PDCO','NIILP','Location','southeast');
lgd = legend('IPRM', 'QOCO','Location','southeast');
set(lgd,'FontName','Times New Roman','FontSize',20,'FontWeight','normal')
xlabel({'log(\tau)'},'FontSize',20,'FontWeight','bold');
ylabel({'\rho_{s}(\tau)'},'FontSize',20,'FontWeight','bold');
saveas(gcf, 'time', 'png');
%title('time');

fclose(fileID1);

