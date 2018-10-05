close all

u= linspace(1,1,10);

n = dirac(u)

stem(f(n))
xlim([-20 20]);
hold on
stem(f(-n))
stem(f(-n+1))
legend('n','-n','-n+1');