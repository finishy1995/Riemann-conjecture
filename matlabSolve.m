%
%  作者：王元恺  日期：2016-10-22
%  Riemann问题matlab处理
%

%%  对比熵修正对Roe格式的影响
% Roe单图
[value1,Roe1]=getPlot('Roe1');
outPlot(value1,Roe1(1,:),'Roe1 速度分布图','Velocity');
print('-dpng','Roe1Velocity');
hold off;
outPlot(value1,Roe1(2,:),'Roe1 密度分布图','Density');
print('-dpng','Roe1Density');
hold off;
outPlot(value1,Roe1(3,:),'Roe1 压力分布图','Pressure');
print('-dpng','Roe1Pressure');
hold off;

% Roe四图比较
[value2,Roe2]=getPlot('Roe2');
[value3,Roe3]=getPlot('Roe3');
[value4,Roe4]=getPlot('Roe4');
[value5,Roe5]=getPlot('Roe5');
figure(1)
subplot(2,2,1);
outPlot(value2,Roe2(1,:),'Roe2 速度分布图','Velocity');
subplot(2,2,2);
outPlot(value3,Roe3(1,:),'Roe3 速度分布图','Velocity');
subplot(2,2,3);
outPlot(value4,Roe4(1,:),'Roe4 速度分布图','Velocity');
subplot(2,2,4);
outPlot(value5,Roe5(1,:),'Roe5 速度分布图','Velocity');
print('-dpng','Roe2345Velocity');
close(1);
figure(2)
subplot(2,2,1);
outPlot(value2,Roe2(2,:),'Roe2 密度分布图','Density');
subplot(2,2,2);
outPlot(value3,Roe3(2,:),'Roe3 密度分布图','Density');
subplot(2,2,3);
outPlot(value4,Roe4(2,:),'Roe4 密度分布图','Density');
subplot(2,2,4);
outPlot(value5,Roe5(2,:),'Roe5 密度分布图','Density');
print('-dpng','Roe2345Density');
close(2);
figure(3)
subplot(2,2,1);
outPlot(value2,Roe2(3,:),'Roe2 压力分布图','Pressure');
subplot(2,2,2);
outPlot(value3,Roe3(3,:),'Roe3 压力分布图','Pressure');
subplot(2,2,3);
outPlot(value4,Roe4(3,:),'Roe4 压力分布图','Pressure');
subplot(2,2,4);
outPlot(value5,Roe5(3,:),'Roe5 压力分布图','Pressure');
print('-dpng','Roe2345Pressure');
close(3);

%%  对比StegerWarming、AUSM和Roe
% SW与Roe对比
[valueS,SW]=getPlot('StegerWarming');
figure(4)
hold on;
plot(value1(1):value1(2):value1(3),Roe1(1,:),'ro','Markersize',4);
plot(valueS(1):valueS(2):valueS(3),SW(1,:),'b','LineWidth',1.5);
title('速度分布图','FontSize',12,'FontWeight','bold','FontName','phong');
xlabel('Position','FontWeight','bold');
ylabel('Velocity','FontWeight','bold');
legend('Roe','StegerWarming');
ylim([0 1]);
print('-dpng','RoeSWVelocity');
hold off;
close(4);
figure(5)
hold on;
plot(value1(1):value1(2):value1(3),Roe1(2,:),'ro','Markersize',4);
plot(valueS(1):valueS(2):valueS(3),SW(2,:),'b','LineWidth',1.5);
title('密度分布图','FontSize',12,'FontWeight','bold','FontName','phong');
xlabel('Position','FontWeight','bold');
ylabel('Density','FontWeight','bold');
legend('Roe','StegerWarming');
ylim([0 1]);
print('-dpng','RoeSWDensity');
hold off;
close(5);
figure(6)
hold on;
plot(value1(1):value1(2):value1(3),Roe1(3,:),'ro','Markersize',4);
plot(valueS(1):valueS(2):valueS(3),SW(3,:),'b','LineWidth',1.5);
title('压力分布图','FontSize',12,'FontWeight','bold','FontName','phong');
xlabel('Position','FontWeight','bold');
ylabel('Pressure','FontWeight','bold');
legend('Roe','StegerWarming');
ylim([0 1]);
print('-dpng','RoeSWPressure');
hold off;
close(6);

% AUSM与Roe对比
[valueA,AUSM]=getPlot('AUSM');
figure(7)
hold on;
plot(value1(1):value1(2):value1(3),Roe1(1,:),'ro','Markersize',4);
plot(valueA(1):valueA(2):valueA(3),AUSM(1,:),'b','LineWidth',1.5);
title('速度分布图','FontSize',12,'FontWeight','bold','FontName','phong');
xlabel('Position','FontWeight','bold');
ylabel('Velocity','FontWeight','bold');
legend('Roe','AUSM');
ylim([0 1]);
print('-dpng','RoeAUSMVelocity');
hold off;
close(7);
figure(8)
hold on;
plot(value1(1):value1(2):value1(3),Roe1(2,:),'ro','Markersize',4);
plot(valueA(1):valueA(2):valueA(3),AUSM(2,:),'b','LineWidth',1.5);
title('密度分布图','FontSize',12,'FontWeight','bold','FontName','phong');
xlabel('Position','FontWeight','bold');
ylabel('Density','FontWeight','bold');
legend('Roe','AUSM');
ylim([0 1]);
print('-dpng','RoeAUSMDensity');
hold off;
close(8);
figure(9)
hold on;
plot(value1(1):value1(2):value1(3),Roe1(3,:),'ro','Markersize',4);
plot(valueA(1):valueA(2):valueA(3),AUSM(3,:),'b','LineWidth',1.5);
title('压力分布图','FontSize',12,'FontWeight','bold','FontName','phong');
xlabel('Position','FontWeight','bold');
ylabel('Pressure','FontWeight','bold');
legend('Roe','AUSM');
ylim([0 1]);
print('-dpng','RoeAUSMPressure');
hold off;
close(9);
