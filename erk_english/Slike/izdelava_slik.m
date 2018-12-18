fontsize = 40;
linewidth = 3;
%% slika 1

k = 0:0.1:10;
h2_amp = 1 .*(k - 1)./(k + 1);

figure(1);
set(gcf, 'Position', get(0, 'Screensize'));
plot(k,h2_amp, 'LineWidth',linewidth, 'Color', [215,25,28]/255);
legend('C_{2}', 'Location','northwest');
set(gca,'FontSize',fontsize);
grid on;
xlabel('$k$','interpreter', 'latex','FontSize',fontsize);
ylabel('$\varepsilon / rad$','interpreter', 'latex','FontSize',fontsize);
saveas(gcf,'amp','epsc')
clear k h2_amp
%% slika 2
fis = 0:0.1:90;
C0 = fis./2;
C2 = 180/pi .* tand(fis./2);
fi2 = fis;
figure(2);
set(gcf, 'Position', get(0, 'Screensize'));
plot(fis,C0,fis,C2,fis,fi2,'Linewidth',linewidth);
xlim([0,90]);
grid on;
set(gca,'FontSize',fontsize);
xlabel('$\mathrm{\varphi_{s}/^\circ}$','Interpreter','latex','FontSize',fontsize);
ylabel('$\varepsilon/^\circ,\varphi/^\circ$','Interpreter','latex','FontSize',fontsize);
leg1 = legend('$C_0$', '$C_2$','$\varphi_2$');
set(leg1,'Interpreter','latex');
set(leg1,'location','NorthWest');
saveas(gcf,'fis','epsc')

clear C0 C2 fi2 fis leg1
%% slika 3

off = -5:0.1:5;
h2off(off<-1) =  (2 +off(off<-1).^-1);
h2off(abs(off)<=1)= -off(abs(off)<=1);
h2off(off>1) = -(2 -off(off>1).^-1);
figure(3);
axes1 = axes('Parent',gcf);
set(gcf, 'Position', get(0, 'Screensize'));
plot(off,h2off,'-r', 'linewidth',linewidth);
grid on;
legend('C_1');
set(axes1,'FontSize',40,'XTick',[-5 -4 -3 -2 -1 0 1 2 3 4 5],'XTickLabel',...
    {'-5','-4','-3','-2','-1','0','1','2','3','4','5'});
set(gca,'FontSize',fontsize);
xlabel('A_0 / A_1','FontSize',fontsize);
ylabel('$\varepsilon / rad$','interpreter', 'latex','FontSize',fontsize);
saveas(gcf,'off','epsc')
clear off h2off

%% slika 4
theta = linspace(0,360-360/4096,4096);
dc = 0:0.05:2;

t1 = zeros(length(dc),length(theta));
for i = 1:length(dc)
    t1(i,:)=theta;
end
Sin = sind(t1)+dc'*cosd(theta);
Cos = cosd(t1)+dc'*cosd(theta);
err = atan2d(Sin,Cos) - atan2d(sind(t1),cosd(t1));
err(err < -180) = err(err < -180) +360;
err(err >  180) = err(err >  180) -360;

L=length(err);
n=L;
    Y=fft(err,n,2);
    P2=abs(Y/n);
    phs=angle(fftshift(Y,2));
    P1=P2(:,1:n/2+1);
    P1(:,2:end-1)=2.*P1(:,2:end-1);
    phase=[phs(:,n/2+1:end),phs(:,1)];
fftErr{1}=P1;
fftErr{2}=phase*180/pi;
c2 = fftErr{1}(:,3);
a2 = c2.*cosd(fftErr{2}(:,3));
b2 = c2.*sind(fftErr{2}(:,3));
c0 = fftErr{1}(:,1);

figure(4)
set(gcf, 'Position', get(0, 'Screensize'));
plot(dc,c0, 'LineWidth',linewidth, 'Color', [215,25,28]/255);
hold on
plot(dc, a2, 'LineWidth',linewidth, 'Color', [253,174,97]/255);
plot(dc, b2, 'LineWidth',linewidth, 'Color', [44,123,182]/255);

legend('C_{0}','C_{2c}', 'C_{2s}', 'Location','northwest')
grid on
set(gca,'FontSize',fontsize);
xlabel('$\Delta_c/A_1$','interpreter', 'latex','FontSize',fontsize)
ylabel('$\varepsilon / ^\circ$','interpreter', 'latex','FontSize',fontsize)
saveas(gcf,'dc','epsc')
clear Cos dc err Fs i L n Sin T t1 theta a2 b2 c0 c2 fftErr P1 P2 phase phs Y
%% slika 5
theta = linspace(0,360-360/4096,4096);
k = 1.1;
error = atan2d(k.*sind(theta), cosd(theta))-atan2d(sind(theta), cosd(theta));
error(error < -180) = error(error < -180) +360;
error(error >  180) = error(error >  180) -360;
predict_e = 0.*theta;

for n = 1:15
    predict_e = predict_e +180/pi .* 1/n .* ((k-1)/(k+1))^n .*sind(2*n.*theta);
end

figure(5)
set(gcf, 'Position', get(0, 'Screensize'));
plot(theta, predict_e-error,'LineWidth',linewidth)
xlim([0,360])
grid on
set(gca,'FontSize',fontsize);
xlabel('$\theta/ ^\circ$','interpreter', 'latex','FontSize',fontsize)
ylabel('$/ ^\circ$','interpreter', 'latex','FontSize',fontsize)
saveas(gcf,'razlika_amp','epsc')    
clear error k n predict_e theta fontsize

