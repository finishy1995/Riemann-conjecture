%
%  ���ߣ���Ԫ��  ���ڣ�2016-10-20
%  ���ͼ��
%

%%  ��ȡͼ������
fid=fopen('Roe.txt','r');
value = [];
tline = fgetl(fid);
tline = str2num(tline);
value = tline;
Roe = [];
for i = 1:3
    tline = fgetl(fid);
    tline = str2num(tline);
    Roe = [Roe;tline];
end
fclose(fid);

%%  ���ͼ��
outPlot(value,Roe(1,:),'RoeVelocity','Roe �ٶȷֲ�ͼ','Velocity');
outPlot(value,Roe(2,:),'RoeDensity','Roe �ܶȷֲ�ͼ','Density');
outPlot(value,Roe(3,:),'RoePressure','Roe ѹ���ֲ�ͼ','Pressure');
