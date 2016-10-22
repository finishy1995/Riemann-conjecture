%
%  作者：王元恺  日期：2016-10-20
%  输出图像
%

%%  获取图像数据
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

%%  输出图像
outPlot(value,Roe(1,:),'RoeVelocity','Roe 速度分布图','Velocity');
outPlot(value,Roe(2,:),'RoeDensity','Roe 密度分布图','Density');
outPlot(value,Roe(3,:),'RoePressure','Roe 压力分布图','Pressure');
