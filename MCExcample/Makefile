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

FASTJET=`$(FASTJETINSTALLPATH)fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`

CHECK_FULLJRDIR=printf "\nFULLJRDIR is set to: \'$(FULLJRDIR)\'\n   !!!IF FULLJRDIR IS EMPTY THE BUILD WILL FAIL!!!\n   Run 'source setStudyMultEnv.sh' in top level StudyMult directory to fix.\n\n"

CHECK_FASTJETINSTALLPATH=printf "\nFASTJETINSTALLPATH is set to: \'$(FASTJETINSTALLPATH)\'\n   !!!IF FASTJETINSTALLPATH IS EMPTY THE BUILD WILL FAIL!!!\n   Run 'source setStudyMultEnv.sh' in top level StudyMult directory to fix.\n\n"

CHECK_ROOUNFPATH=printf "\nROOUNFPATH is set to: \'$(ROOUNFPATH)\'\n   !!!IF ROOUNFPATH IS EMPTY THE BUILD WILL FAIL!!!\n   Run 'source setStudyMultEnv.sh' in top level StudyMult directory to fix.\n\n"

all: mkdirBin mkdirPdf mkdirOutput checkFullJRDir checkFastJetInstallPath checkRooUnfPath bin/processAndRecoPYT8.exe bin/makeRawRAASpectra.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdf:
	$(MKDIR_PDF)

checkFullJRDir:
	@$(CHECK_FULLJRDIR)

checkFastJetInstallPath:
	@$(CHECK_FASTJETINSTALLPATH)

checkRooUnfPath:
	@$(CHECK_ROOUNFPATH)

bin/processAndRecoPYT8.exe: src/processAndRecoPYT8.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(FASTJET) $(INCLUDE)  -o bin/processAndRecoPYT8.exe src/processAndRecoPYT8.C

bin/makeRawRAASpectra.exe: src/makeRawRAASpectra.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(FASTJET) $(INCLUDE) -I $(ROOUNFPATH) -L $(ROOUNFPATH) -lRooUnfold -o bin/makeRawRAASpectra.exe src/makeRawRAASpectra.C

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