%
%  ���ߣ���Ԫ��  ���ڣ�2016-10-22
%  �����ļ�����
%

function [value,Roe]=getPlot(filename)
    %%  ��ȡͼ������
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
