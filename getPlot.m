%
%  ���ߣ���Ԫ��  ���ڣ�2016-10-22
%  �����ļ����ݲ����ͼ��
%

function getPlot(filename,bool)
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

    %%  ���ͼ��
    outPlot(value,Roe(1,:),[filename,'Velocity'],'Roe �ٶȷֲ�ͼ','Velocity',bool);
    outPlot(value,Roe(2,:),[filename,'Density'],'Roe �ܶȷֲ�ͼ','Density',bool);
    outPlot(value,Roe(3,:),[filename,'Pressure'],'Roe ѹ���ֲ�ͼ','Pressure',bool);
end
