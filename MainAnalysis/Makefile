CXX = g++
#O3 for max optimization (go to 0 for debug)
CXXFLAGS = -Wall -Werror -O0 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g
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

all: mkdirBin mkdirPdf mkdirOutput checkFullJRDir bin/makeRawRAATree.exe bin/makeJetResponseTree.exe bin/makeDeriveFlatResponse.exe bin/validateJetResponse.exe bin/plotJetResponse.exe bin/makeMuonFakeCheck.exe bin/checkGoodJetBadJetPF.exe bin/processRawData.exe bin/unfoldRawData.exe bin/deriveDijetWeights.exe bin/checkGenHydSpectra.exe bin/plotUnfoldedSpectra.exe bin/plotUnfoldedAll.exe bin/checkUnfoldTermination.exe bin/combineFiles.exe bin/testRooUnfResUnderflow.exe bin/checkFileNClassContents.exe #bin/plotUnfoldedRAA.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdf:
	$(MKDIR_PDF)

checkFullJRDir:
	@$(CHECK_FULLJRDIR)

bin/makeRawRAATree.exe: src/makeRawRAATree.C
	$(CXX) $(CXXFLAGS)  src/makeRawRAATree.C  $(ROOT) $(INCLUDE) -o bin/makeRawRAATree.exe

bin/makeJetResponseTree.exe: src/makeJetResponseTree.C
	$(CXX) $(CXXFLAGS) src/makeJetResponseTree.C  $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/makeJetResponseTree.exe

bin/makeDeriveFlatResponse.exe: src/makeDeriveFlatResponse.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(INCLUDE) -o bin/makeDeriveFlatResponse.exe src/makeDeriveFlatResponse.C

bin/validateJetResponse.exe: src/validateJetResponse.C
	$(CXX) $(CXXFLAGS) src/validateJetResponse.C $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/validateJetResponse.exe

bin/plotJetResponse.exe: src/plotJetResponse.C
	$(CXX) $(CXXFLAGS) src/plotJetResponse.C $(ROOT) $(INCLUDE) -o bin/plotJetResponse.exe

bin/plotUnfoldedSpectra.exe: src/plotUnfoldedSpectra.C
	$(CXX) $(CXXFLAGS) src/plotUnfoldedSpectra.C $(ROOT) $(INCLUDE) -o bin/plotUnfoldedSpectra.exe

bin/plotUnfoldedAll.exe: src/plotUnfoldedAll.C
	$(CXX) $(CXXFLAGS) src/plotUnfoldedAll.C $(ROOT) $(INCLUDE) -o bin/plotUnfoldedAll.exe

bin/plotUnfoldedRAA.exe: src/plotUnfoldedRAA.C
	$(CXX) $(CXXFLAGS) src/plotUnfoldedRAA.C $(ROOT) $(INCLUDE) -o bin/plotUnfoldedRAA.exe

bin/makeMuonFakeCheck.exe: src/makeMuonFakeCheck.C
	$(CXX) $(CXXFLAGS) src/makeMuonFakeCheck.C $(ROOT) $(INCLUDE) -o bin/makeMuonFakeCheck.exe

bin/checkGoodJetBadJetPF.exe: src/checkGoodJetBadJetPF.C
	$(CXX) $(CXXFLAGS) src/checkGoodJetBadJetPF.C $(ROOT) $(INCLUDE) -o bin/checkGoodJetBadJetPF.exe

bin/checkGenHydSpectra.exe: src/checkGenHydSpectra.C
	$(CXX) $(CXXFLAGS) src/checkGenHydSpectra.C $(ROOT) $(INCLUDE) -o bin/checkGenHydSpectra.exe

bin/processRawData.exe: src/processRawData.C
	$(CXX) $(CXXFLAGS) src/processRawData.C $(ROOT) $(INCLUDE) -o bin/processRawData.exe

bin/unfoldRawData.exe: src/unfoldRawData.C
	$(CXX) $(CXXFLAGS) src/unfoldRawData.C $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/unfoldRawData.exe

bin/deriveDijetWeights.exe: src/deriveDijetWeights.C
	$(CXX) $(CXXFLAGS) src/deriveDijetWeights.C $(ROOT) $(INCLUDE) -o bin/deriveDijetWeights.exe

bin/checkUnfoldTermination.exe: src/checkUnfoldTermination.C
	$(CXX) $(CXXFLAGS) src/checkUnfoldTermination.C $(ROOT) $(INCLUDE) -o bin/checkUnfoldTermination.exe

bin/combineFiles.exe: src/combineFiles.C
	$(CXX) $(CXXFLAGS) src/combineFiles.C $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/combineFiles.exe

bin/testRooUnfResUnderflow.exe: src/testRooUnfResUnderflow.C
	$(CXX) $(CXXFLAGS) src/testRooUnfResUnderflow.C $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/testRooUnfResUnderflow.exe

bin/checkFileNClassContents.exe: src/checkFileNClassContents.C
	$(CXX) $(CXXFLAGS) src/checkFileNClassContents.C $(ROOT) $(INCLUDE) -o bin/checkFileNClassContents.exe

clean:
	rm -f ./*~
	rm -f ./#*#
	rm -f bash/*~
	rm -f bash/#*#
	rm -f bin/*.exe
	rm -rf bin
	rm -f configs/*~
	rm -f configs/#*#
	rm -f include/*~
	rm -f include/#*#
	rm -f src/*~
	rm -f src/#*#
