%
%  ���ߣ���Ԫ��  ���ڣ�2016-10-22
%  ���ͼ��
%

function outPlot(value,matrix,filename,plotname,yname,bool)
    plot(value(1):value(2):value(3),matrix);
    if bool==1
        title([plotname,'  Tol=',num2str(value(4))],'FontSize',12,'FontWeight','bold','FontName','phong');
    else
        title(plotname,'FontSize',12,'FontWeight','bold','FontName','phong');
    end
    xlabel('Position','FontWeight','bold');
    ylabel(yname,'FontWeight','bold');
    print('-dpng',filename);
end
