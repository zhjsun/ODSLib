# 目录设置
SrcDir = ./src
BuildDir = ./build
# 中间文件
Objects =  	$(BuildDir)/MainRun.o
# 目标文件
Target = runner
# 编译器设置
CC = clang++
# 编译选项设置
Compiler = $(Debug) -std=c++11
Debug = -g
Release = -O3
# OpenMP并行设置
OpenMP = -fopenmp
# 包含头文件目录
Include = -I .
# 连接产生目标文件
$(Target) : $(Objects)
	$(CC) -o $(Target) $(Objects)
# 各源文件编译
$(BuildDir)/MainRun.o : $(SrcDir)/MainRun.cpp
	$(CC) -c $(SrcDir)/MainRun.cpp -o $(BuildDir)/MainRun.o $(Include) $(Compiler)
# 清除中间文件
clean : 
	rm  $(Objects) $(Target)
rebuild :
	make clean
	make
