$(info INC: $(INC))  
$(info PYTHIAINC: $(PYTHIAINC))  
$(info PYTHIALIB: $(PYTHIALIB))  
$(info ROOTLIBS: $(ROOTLIBS))

# 设置编译器
CXX = g++

# 设置 ROOT 头文件路径和库路径
#ROOTSYS = /home/mihu/home/root  # 替换为 ROOT 的实际安装路径
INC += -I$(ROOTSYS)/include 
LIB += -L$(ROOTSYS)/lib    

# 设置需要链接的 ROOT 库
ROOTLIBS = $(shell root-config --libs)
ROOTCFLAGS = $(shell root-config --cflags)

# 设置 PYTHIA8 头文件路径和库路径
#PYTHIA8 = /home/mihu/home/pythia8311
PYTHIAINC = -I$(PYTHIA8)/include
PYTHIALIB = -L$(PYTHIA8)/lib

# 设置编译选项
CXXFLAGS = -g $(INC) $(PYTHIAINC) $(ROOTCFLAGS)

# 设置要编译的源文件
SRC = hello.c

# 设置输出的可执行文件
EXE = hello.exe

# 默认目标
all: $(EXE)

# 链接目标
$(EXE): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(EXE) $(PYTHIALIB) -lpythia8 $(ROOTLIBS)

# 清理目标
clean:
	rm -f $(EXE) *.o *~

# 伪目标
.PHONY: all clean
