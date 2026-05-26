
#
# Makefile
#

CC=nvcc
EXE=.\main.exe
EXT=cu
OBJ=obj

# gigaByte PC

OPTIONS= -ccbin "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\VC\Tools\MSVC\14.29.30133\bin\Hostx86\x64"

#Aorus PC
#OPTIONS= -ccbin "C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Tools\MSVC\14.29.30133\bin\Hostx86\x64"

# naoki PC
#OPTIONS= -ccbin "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Tools\MSVC\14.44.35207\bin\Hostx86\x64" -Xcompiler "/wd 4819" -diag-suppress 221 --std c++17

# MSI PC
#OPTIONS= -ccbin "C:\Program Files\Microsoft Visual Studio\18\Enterprise\VC\Tools\MSVC\14.50.35717\bin\Hostx86\x64" -Xcompiler /Zc:preprocessor -Xcompiler "/wd 4819" --std c++17 -diag-suppress 221

# # for Ubuntu
# # 1. Do first "make cp"
# CC=g++
# OPTIONS= -D UBUNTU
# EXE=./main
# EXT=cpp
# OBJ=o

SRCS=main.${EXT} PCSData.${EXT} LambertWhara.${EXT} adaptiveModel.${EXT} substring.${EXT} stringModel.${EXT} stringGroup.${EXT} globalConstants.${EXT} searchHistoryAM.${EXT}
OBJS=$(SRCS:.${EXT}=.${OBJ})
DEPS=$(SRCS:.${EXT}=.d)

all: ${OBJS}
	${CC} ${OPTIONS} -o main ${OBJS}
%.${OBJ}: %.${EXT}
	${CC} ${OPTIONS} -c -o $@ $<
%.d: %.${EXT}
	${CC} -M ${OPTIONS} $< > $@

# for Ubuntu
cp:
	chmod +x cp2cpp.sh
	./cp2cpp.sh

go:
	${EXE}

clean:
	powershell rm *.exe
	powershell rm *.obj
	powershell rm *.lib
	powershell rm *.exp
	powershell rm *.d
	powershell rm *~
	powershell rm Graphs/*~
	powershell rm Dat/*~
	powershell rm Doc/*~
	powershell rm Python/*~

cleanUbuntu:
	rm -rf *~

# MSI PC, Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope Process
zip:
	powershell Compress-Archive -Update -Path ".\*.cu", "*.hh", ".\Makefile", "*.ps1", "*.odt", "*.cc", "*.h", "*.sh", "Graphs", "Dat", "Doc", "Python"  -DestinationPath ".\stringGroupCpp.zip"

tar:
	tar cvjf stringGroupCpp.tar.bz2 *.cu *.hh Makefile *.ps1 *.odt *.cc *.h *.sh Graphs/ Dat/ Doc/ Python/

ls:
	powershell Set-ExecutionPolicy -Scope CurrentUser Unrestricted
	powershell .\ls.ps1

txt:
	powershell get-content main.cu, PCSData.hh, PCSData.cu, HornerHara.hh, LambertWhara.hh, LambertWhara.cu, globalConstants.hh, adaptiveModel.hh, adaptiveModel.cu, bisection.hh, substring.hh, substring.cu, stringModel.hh, stringModel.cu, stringGroup.hh, stringGroup.cu  > program.txt

dep: ${DEPS}

-include ${DEPS}

