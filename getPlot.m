%
%  作者：王元恺  日期：2016-10-22
%  读入文件数据并输出图像
%

function getPlot(filename,bool)
    %%  获取图像数据
    fid=fopen([filename,'.txt'],'r');
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
    outPlot(value,Roe(1,:),[filename,'Velocity'],'Roe 速度分布图','Velocity',bool);
    outPlot(value,Roe(2,:),[filename,'Density'],'Roe 密度分布图','Density',bool);
    outPlot(value,Roe(3,:),[filename,'Pressure'],'Roe 压力分布图','Pressure',bool);
end
