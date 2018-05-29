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

ROOUNF=-I $(ROOUNFOLDPATH) -L $(ROOUNFOLDPATH) -lRooUnfold

all: mkdirBin mkdirPdf mkdirOutput checkFullJRDir bin/makeRawRAATree.exe bin/makeJetResponseTree.exe bin/validateJetResponse.exe bin/plotJetResponse.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdf:
	$(MKDIR_PDF)

checkFullJRDir:
	@$(CHECK_FULLJRDIR)

bin/makeRawRAATree.exe: src/makeRawRAATree.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(INCLUDE) -o bin/makeRawRAATree.exe src/makeRawRAATree.C

bin/makeJetResponseTree.exe: src/makeJetResponseTree.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/makeJetResponseTree.exe src/makeJetResponseTree.C

bin/validateJetResponse.exe: src/validateJetResponse.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/validateJetResponse.exe src/validateJetResponse.C

bin/plotJetResponse.exe: src/plotJetResponse.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(INCLUDE) -o bin/plotJetResponse.exe src/plotJetResponse.C

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