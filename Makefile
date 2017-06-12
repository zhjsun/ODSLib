# 目录设置
SrcDir = ./src
BuildDir = ./build
# 中间文件
Objects =  	$(BuildDir)/MainRun.o \
			$(BuildDir)/ODS.o
# 目标文件
Target = bin/runner
# 编译器设置
CC = g++
# 编译选项设置
Compiler = $(Release) -std=c++11
Debug = -g -Og
Release = -O3
# OpenMP并行设置
OpenMP = -fopenmp
# 包含头文件目录
Include = -I ./include \
		  -I ./src
# 连接产生目标文件
$(Target) : $(Objects)
	$(CC) -o $(Target) $(Objects)
# 各源文件编译
$(BuildDir)/MainRun.o : $(SrcDir)/MainRun.cpp
	$(CC) -c $(SrcDir)/MainRun.cpp -o $(BuildDir)/MainRun.o $(Include) $(Compiler)
$(BuildDir)/ODS.o : $(SrcDir)/ODSLib/ODS.cpp
	$(CC) -c $(SrcDir)/ODSLib/ODS.cpp -o $(BuildDir)/ODS.o $(Include) $(Compiler)
# 清除中间文件
clean : 
	rm  $(Objects) $(Target)
rebuild :
	make clean
	make
