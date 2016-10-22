function outPlot(value,matrix,filename,plotname,yname)
    plot(value(1):value(2):value(3),matrix);
    title(plotname,'FontSize',12,'FontWeight','bold','FontName','phong');
    xlabel('Position','FontWeight','bold');
    ylabel(yname,'FontWeight','bold');
    print('-dpng',filename);
end