SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .c .o
CC = cc
CFLAGS = -g -DDEBUG -D_CONSOLE -O3 -funroll-loops
LDFLAGS = -lm -lexpat

OBJS = base64.o BN.o BuildMS2DB.o ChargeState.o CMemLeak.o Errors.o ExonGraphAlign.o \
	FreeMod.o IonScoring.o \
	LDA.o main.o Mods.o MS2DB.o ParentMass.o ParseInput.o ParseXML.o PValue.o \
	Run.o Score.o Scorpion.o SNP.o Spectrum.o Spliced.o \
	SpliceDB.o SpliceScan.o SVM.o Tagger.o Trie.o Utils.o 
HDRS = base64.h BN.h BuildMS2DB.h ChargeState.h CMemLeak.h Errors.h ExonGraphAlign.h FreeMod.h \
	Inspect.h IonScoring.h LDA.h Mods.h MS2DB.h ParentMass.h ParseInput.h ParseXML.h PValue.h \
	Run.h Score.h Scorpion.h \
	SNP.h Spectrum.h Spliced.h SpliceDB.h SpliceScan.h SVM.h Tagger.h \
	Trie.h Utils.h 

EXE = inspect

all: $(EXE)

.c.o: $(HDRS)
	$(CC) $(CFLAGS) -c $<

$(EXE): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(EXE) $(OBJS)

clean-objs:
	@rm -f $(OBJS) $(MPI_OBJS)

clean-exe:
	@rm -f $(EXE) $(MPI_EXE) 

clean: clean-objs clean-exe
