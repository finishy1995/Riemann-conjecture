%
%  ���ߣ���Ԫ��  ���ڣ�2016-10-22
%  Riemann����matlab����
%

%%  �Ա���������Roe��ʽ��Ӱ��
% Roe��ͼ
[value1,Roe1]=getPlot('Roe1');
outPlot(value1,Roe1(1,:),'Roe1 �ٶȷֲ�ͼ','Velocity');
print('-dpng','Roe1Velocity');
hold off;
outPlot(value1,Roe1(2,:),'Roe1 �ܶȷֲ�ͼ','Density');
print('-dpng','Roe1Density');
hold off;
outPlot(value1,Roe1(3,:),'Roe1 ѹ���ֲ�ͼ','Pressure');
print('-dpng','Roe1Pressure');
hold off;

% Roe��ͼ�Ƚ�
[value2,Roe2]=getPlot('Roe2');
[value3,Roe3]=getPlot('Roe3');
[value4,Roe4]=getPlot('Roe4');
[value5,Roe5]=getPlot('Roe5');
figure(1)
subplot(2,2,1);
outPlot(value2,Roe2(1,:),'Roe2 �ٶȷֲ�ͼ','Velocity');
subplot(2,2,2);
outPlot(value3,Roe3(1,:),'Roe3 �ٶȷֲ�ͼ','Velocity');
subplot(2,2,3);
outPlot(value4,Roe4(1,:),'Roe4 �ٶȷֲ�ͼ','Velocity');
subplot(2,2,4);
outPlot(value5,Roe5(1,:),'Roe5 �ٶȷֲ�ͼ','Velocity');
print('-dpng','Roe2345Velocity');
close(1);
figure(2)
subplot(2,2,1);
outPlot(value2,Roe2(2,:),'Roe2 �ܶȷֲ�ͼ','Density');
subplot(2,2,2);
outPlot(value3,Roe3(2,:),'Roe3 �ܶȷֲ�ͼ','Density');
subplot(2,2,3);
outPlot(value4,Roe4(2,:),'Roe4 �ܶȷֲ�ͼ','Density');
subplot(2,2,4);
outPlot(value5,Roe5(2,:),'Roe5 �ܶȷֲ�ͼ','Density');
print('-dpng','Roe2345Density');
close(2);
figure(3)
subplot(2,2,1);
outPlot(value2,Roe2(3,:),'Roe2 ѹ���ֲ�ͼ','Pressure');
subplot(2,2,2);
outPlot(value3,Roe3(3,:),'Roe3 ѹ���ֲ�ͼ','Pressure');
subplot(2,2,3);
outPlot(value4,Roe4(3,:),'Roe4 ѹ���ֲ�ͼ','Pressure');
subplot(2,2,4);
outPlot(value5,Roe5(3,:),'Roe5 ѹ���ֲ�ͼ','Pressure');
print('-dpng','Roe2345Pressure');
close(3);

%%  �Ա�StegerWarming��AUSM��Roe
% SW��Roe�Ա�
[valueS,SW]=getPlot('StegerWarming');
figure(4)
hold on;
plot(value1(1):value1(2):value1(3),Roe1(1,:),'ro','Markersize',4);
plot(valueS(1):valueS(2):valueS(3),SW(1,:),'b','LineWidth',1.5);
title('�ٶȷֲ�ͼ','FontSize',12,'FontWeight','bold','FontName','phong');
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
title('�ܶȷֲ�ͼ','FontSize',12,'FontWeight','bold','FontName','phong');
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
title('ѹ���ֲ�ͼ','FontSize',12,'FontWeight','bold','FontName','phong');
xlabel('Position','FontWeight','bold');
ylabel('Pressure','FontWeight','bold');
legend('Roe','StegerWarming');
ylim([0 1]);
print('-dpng','RoeSWPressure');
hold off;
close(6);

% AUSM��Roe�Ա�
[valueA,AUSM]=getPlot('AUSM');
figure(7)
hold on;
plot(value1(1):value1(2):value1(3),Roe1(1,:),'ro','Markersize',4);
plot(valueA(1):valueA(2):valueA(3),AUSM(1,:),'b','LineWidth',1.5);
title('�ٶȷֲ�ͼ','FontSize',12,'FontWeight','bold','FontName','phong');
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
title('�ܶȷֲ�ͼ','FontSize',12,'FontWeight','bold','FontName','phong');
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
title('ѹ���ֲ�ͼ','FontSize',12,'FontWeight','bold','FontName','phong');
xlabel('Position','FontWeight','bold');
ylabel('Pressure','FontWeight','bold');
legend('Roe','AUSM');
ylim([0 1]);
print('-dpng','RoeAUSMPressure');
hold off;
close(9);
