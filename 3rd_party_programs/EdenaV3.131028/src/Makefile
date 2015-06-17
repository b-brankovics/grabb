CCC 	=	g++
CFLAGS	=	-O2 -ansi -pthread

all: edena

node.o: node.cpp
	$(CCC) -o node.o -c $(CFLAGS) node.cpp
	
overlapGraph.o: overlapGraph.cpp
	$(CCC) -o overlapGraph.o -c $(CFLAGS) overlapGraph.cpp
	
readsStorage.o: readsStorage.cpp
	$(CCC) -o readsStorage.o -c $(CFLAGS) readsStorage.cpp
	
stat.o: stat.cpp
	$(CCC) -o stat.o -c $(CFLAGS) stat.cpp
	
NodeIt.o: NodeIt.cpp
	$(CCC) -o NodeIt.o -c $(CFLAGS) NodeIt.cpp
	
Pairing.o: Pairing.cpp
	$(CCC) -o Pairing.o -c $(CFLAGS) Pairing.cpp
	
readsLayout.o: readsLayout.cpp
	$(CCC) -o readsLayout.o -c $(CFLAGS) readsLayout.cpp
	
globalFunc.o: globalFunc.cpp
	$(CCC) -o globalFunc.o -c $(CFLAGS) globalFunc.cpp
	
DevShell.o: DevShell.cpp
	$(CCC) -o DevShell.o -c $(CFLAGS) DevShell.cpp
	
BeamSearchTree.o: BeamSearchTree.cpp
	$(CCC) -o BeamSearchTree.o -c $(CFLAGS) BeamSearchTree.cpp

BackwardsWalker.o: BackwardsWalker.cpp
	$(CCC) -o BackwardsWalker.o -c $(CFLAGS) BackwardsWalker.cpp

PEMatches.o: PEMatches.cpp
	$(CCC) -o PEMatches.o -c $(CFLAGS) PEMatches.cpp
	
ActualDistribution.o: ActualDistribution.cpp
	$(CCC) -o ActualDistribution.o -c $(CFLAGS) ActualDistribution.cpp

logWriter.o: logWriter.cpp
	$(CCC) -o logWriter.o -c $(CFLAGS) logWriter.cpp

Param.o: Param.cpp
	$(CCC) -o Param.o -c $(CFLAGS) Param.cpp
	
edena:  node.o overlapGraph.o readsStorage.o stat.o BackwardsWalker.o \
	NodeIt.o Pairing.o readsLayout.o globalFunc.o DevShell.o \
	BeamSearchTree.o PEMatches.o ActualDistribution.o logWriter.o \
	Param.o
	$(CCC) -o edena $(CFLAGS) main.cpp \
	node.o overlapGraph.o readsStorage.o stat.o BackwardsWalker.o \
	NodeIt.o Pairing.o readsLayout.o globalFunc.o DevShell.o \
	BeamSearchTree.o PEMatches.o ActualDistribution.o logWriter.o \
	Param.o
	
clean:
	rm node.o overlapGraph.o readsStorage.o stat.o BackwardsWalker.o \
	NodeIt.o Pairing.o \
	readsLayout.o globalFunc.o DevShell.o \
	BeamSearchTree.o PEMatches.o ActualDistribution.o logWriter.o \
	Param.o
	
