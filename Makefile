# 中间文件
Objects =  	MainRun.o \
			ODSLib.o
# 编译器设置
CC = g++
# 编译选项设置
Compiler = $(Release) -std=c++11
Debug = -g -Og
Release = -O3
# OpenMP并行设置
OpenMP = -fopenmp
# 包含头文件目录
Include = -I .\
	  -I Eigen
# 连接产生目标文件
Test : $(Objects)
	$(CC) -o Test.exe $(Objects) $(OpenMP)
# 各源文件编译
MainRun.o : MainRun.cpp
	$(CC) -c MainRun.cpp $(Include) $(Compiler)
ODSLib.o : ODSLib.cpp
	$(CC) -c ODSLib.cpp $(Include) $(Compiler)
# 清除中间文件
clean : 
	rm  $(Objects) Test.exe
rebuild :
	make clean
	make
