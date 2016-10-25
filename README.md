# Riemann-conjecture
Number format to solve Riemann conjecture(C++ & Matlab)

author: David Wang  王元恺

ShangHai Jiao Tong University  上海交通大学

先编译运行main.cpp，后运行matlabSolve.m
--------------------------------------

c++执行文件 main.cpp

抽象基类 NumFormat.hpp 定义方法的初始值和公共方法

普通类 Roe.hpp 继承自NumFormat.hpp

普通类 StegerWarming.hpp 继承自NumFormat.hpp

普通类 AUSM.hpp 继承自NumFormat.hpp

常用方法：

    {objectname}.setFilePath(FILEPATH); 设置输出路径

    {objectname}.solve(); 执行运算方法

    {objectname}.output("Roe1.txt"); 输出至txt文件

matlab执行文件 matlabSolve.m
