A = [];
b = [];
Aeq = [];
beq = [];
lb = 1e-8;
ub = 1e-4;
nonlcon = [];
options = optimoptions('fmincon','Display','iter');
D0 = 5e-6;
D = fmincon(@costfun,D0,A,b,Aeq,beq,lb,ub,nonlcon,options);
syms n t
L= 0.04;
f_1(t)=(8/((2*n+1)^2*pi^2))*exp((-D*(2*n+1)^2*pi^2*t)/(4*L^2));
M_1(t)=1-symsum(f_1(t),n,0,Inf);
T=exp([0:0.2:ceil(log(7200))]);
mt_minf_1=double(M_1(T));
figure()
xlabel('Time (hours)','Interpreter','latex')
ylabel('$$M_{t}/M_{\infty} (\%)$$','Interpreter','latex');
set (gca,'FontSize', 16)
hold on
p = plot(Time,100*Diffusion_Fraction,"ko ");
p.MarkerFaceColor = [0 0 0];
p.MarkerSize = 8;
plot(T/3600,100*mt_minf_1,'r','LineWidth',1)
legend('Trypsin Inhibitor', 'Fitted Data','Interpreter','latex');
text_show = ['$$D = ' num2str(round(D/1e-6,2)) ' \times 10^{-6} cm^{2}/s$$'];
text(1,30,text_show, 'Interpreter','latex', 'FontSize', 16)
box on;
hold off
%%%%%
function cost = costfun(D)
global Time Diffusion_Fraction
syms n t
L= 0.04;
f_1(t)=(8/((2*n+1)^2*pi^2))*exp((-D*(2*n+1)^2*pi^2*t)/(4*L^2));
M_1(t)=1-symsum(f_1(t),n,0,Inf);
T=3600*Time;
mt_minf_1=double(M_1(T));
output = mt_minf_1';
cost = 1e4*sum((Diffusion_Fraction' - output).^2);
end
