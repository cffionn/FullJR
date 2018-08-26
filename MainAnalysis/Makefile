CXX = g++
CXXFLAGS = -Wall -Werror -O2 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

ifndef FULLJRDIR
$(error FULLJRDIR is not set at all. Do 'source setFullJREnv.sh')
endif

INCLUDE=-I $(FULLJRDIR)
ROOT=`root-config --cflags --glibs`

MKDIR_BIN=mkdir -p $(PWD)/bin
MKDIR_OUTPUT=mkdir -p $(PWD)/output
MKDIR_PDF=mkdir -p $(PWD)/pdfDir

CHECK_FULLJRDIR=printf "\nFULLJRDIR is set to: \'$(FULLJRDIR)\'\n   !!!IF FULLJRDIR IS EMPTY THE BUILD WILL FAIL!!!\n   Run 'source \
setStudyMultEnv.sh' in top level StudyMult directory to fix.\n\n"

ROOUNF=-I $(ROOUNFPATH) -L $(ROOUNFPATH) -lRooUnfold

all: mkdirBin mkdirPdf mkdirOutput checkFullJRDir bin/makeRawRAATree.exe bin/makeJetResponseTree.exe bin/makeDeriveFlatResponse.exe bin/validateJetResponse.exe bin/plotJetResponse.exe bin/makeMuonFakeCheck.exe bin/checkGoodJetBadJetPF.exe bin/processRawData.exe bin/unfoldRawData.exe bin/deriveDijetWeights.exe bin/checkGenHydSpectra.exe bin/plotUnfoldedSpectra.exe bin/checkUnfoldTermination.exe bin/combineFiles.exe bin/testRooUnfResUnderflow.exe bin/checkFileNClassContents.exe

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

bin/makeDeriveFlatResponse.exe: src/makeDeriveFlatResponse.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(INCLUDE) -o bin/makeDeriveFlatResponse.exe src/makeDeriveFlatResponse.C

bin/validateJetResponse.exe: src/validateJetResponse.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/validateJetResponse.exe src/validateJetResponse.C

bin/plotJetResponse.exe: src/plotJetResponse.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(INCLUDE) -o bin/plotJetResponse.exe src/plotJetResponse.C

bin/plotUnfoldedSpectra.exe: src/plotUnfoldedSpectra.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(INCLUDE) -o bin/plotUnfoldedSpectra.exe src/plotUnfoldedSpectra.C

bin/makeMuonFakeCheck.exe: src/makeMuonFakeCheck.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(INCLUDE) -o bin/makeMuonFakeCheck.exe src/makeMuonFakeCheck.C

bin/checkGoodJetBadJetPF.exe: src/checkGoodJetBadJetPF.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(INCLUDE) -o bin/checkGoodJetBadJetPF.exe src/checkGoodJetBadJetPF.C

bin/checkGenHydSpectra.exe: src/checkGenHydSpectra.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(INCLUDE) -o bin/checkGenHydSpectra.exe src/checkGenHydSpectra.C

bin/processRawData.exe: src/processRawData.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(INCLUDE) -o bin/processRawData.exe src/processRawData.C

bin/unfoldRawData.exe: src/unfoldRawData.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/unfoldRawData.exe src/unfoldRawData.C

bin/deriveDijetWeights.exe: src/deriveDijetWeights.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(INCLUDE) -o bin/deriveDijetWeights.exe src/deriveDijetWeights.C

bin/checkUnfoldTermination.exe: src/checkUnfoldTermination.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(INCLUDE) -o bin/checkUnfoldTermination.exe src/checkUnfoldTermination.C

bin/combineFiles.exe: src/combineFiles.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/combineFiles.exe src/combineFiles.C

bin/testRooUnfResUnderflow.exe: src/testRooUnfResUnderflow.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/testRooUnfResUnderflow.exe src/testRooUnfResUnderflow.C

bin/checkFileNClassContents.exe: src/checkFileNClassContents.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(INCLUDE) -o bin/checkFileNClassContents.exe src/checkFileNClassContents.C

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