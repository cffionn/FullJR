CXX = g++
CXXFLAGS = -Wall -Werror -O2 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

INCLUDE=-I $(FULLJRDIR)
ROOT=`root-config --cflags --glibs`

MKDIR_BIN=mkdir -p $(PWD)/bin
MKDIR_OUTPUT=mkdir -p $(PWD)/output
MKDIR_PDF=mkdir -p $(PWD)/pdfDir

CHECK_FULLJRDIR=printf "\nFULLJRDIR is set to: \'$(FULLJRDIR)\'\n   !!!IF FULLJRDIR IS EMPTY THE BUILD WILL FAIL!!!\n   Run 'source \
setStudyMultEnv.sh' in top level StudyMult directory to fix.\n\n"

all: mkdirBin mkdirPdf mkdirOutput checkFullJRDir bin/recreateV2V3Tree.exe bin/recreateV2V3TreeHist.exe bin/plotRecreateV2V3Tree.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdf:
	$(MKDIR_PDF)

checkFullJRDir:
	@$(CHECK_FULLJRDIR)

bin/recreateV2V3Tree.exe: src/recreateV2V3Tree.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(INCLUDE) -o bin/recreateV2V3Tree.exe src/recreateV2V3Tree.C

bin/recreateV2V3TreeHist.exe: src/recreateV2V3TreeHist.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(INCLUDE) -o bin/recreateV2V3TreeHist.exe src/recreateV2V3TreeHist.C

bin/plotRecreateV2V3Tree.exe: src/plotRecreateV2V3Tree.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(INCLUDE) -o bin/plotRecreateV2V3Tree.exe src/plotRecreateV2V3Tree.C


clean:
	rm -f ./*~
	rm -f src/*~
	rm -f include/*~
	rm -f configs/*~
	rm -f ./#*#
	rm -f src/#*#
	rm -f include/#*#
	rm -f configs/#*#
	rm -f bin/*.exe
	rm -rf bin