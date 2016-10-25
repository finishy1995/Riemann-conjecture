%
%  作者：王元恺  日期：2016-10-22
%  读入文件数据
%

function [value,Roe]=getPlot(filename)
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
end
