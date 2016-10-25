%
%  作者：王元恺  日期：2016-10-22
%  输出图像
%

function outPlot(value,matrix,plotname,yname)
    plot(value(1):value(2):value(3),matrix);
    hold on;
    title([plotname,'  Tol=',num2str(value(4))],'FontSize',12,'FontWeight','bold','FontName','phong');
    xlabel('Position','FontWeight','bold');
    ylabel(yname,'FontWeight','bold');
    ylim([0 1]);
end
