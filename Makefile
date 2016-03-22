
ROOTCFLAGS    = $(shell root-config --cflags)

MYCF=""

ifneq ($(findstring m32,$(ROOTCFLAGS)),)
        MYCF= CXXFLAGS=-m32 CFLAGS=-m32
endif

ROOTLIB    = $(shell root-config --glibs) -lMinuit -lTreePlayer -lMathCore

# Get the RestFrames specific compile flags
#RFCFLAGS   = $(shell restframes-config --cxxflags)
#RFLIBS     = $(shell restframes-config --libs)


CXXFLAGS   = -g -I.
CXXFLAGS  += -Wno-long-long -fPIC
CXXFLAGS  += $(shell root-config --cflags)
#CXXFLAGS  += $(RFCFLAGS)

LDFLAGS    =
LDFLAGS   += $(ROOTLIB)
LDFlAGS   += $(RFLIBS)

OBJS       = Root/CollectionTree.o
OBJS      += Root/ChainHelper.o
OBJS      += Root/string_utils.o

OBJS_READ += $(OBJS)
OBJS_READ += util/runStopDiag0L.o

# can use RootCore, but keep like this to make it portable
#OBJS_READ += TopDataPreparation/Root/SampleXsection.o
# CXXFLAGS += -I../RootCoreBin/include
CXXFLAGS += -I ./STop_Ntuple

#  to include the HistoList tool
# CXXFLAGS += -I../TuDoBase
# LDFLAGS += -L../TuDoBase/lib/ -lTuDoBase

%.o: %.cxx
	g++ -c $(CXXFLAGS) -o $@ $< 

all: runStopDiag0L

runStopDiag0L: $(OBJS_READ)
	g++ $(CXXFLAGS) -o runStopDiag0L $(OBJS_READ) $(LDFLAGS)

clean:
	rm -rf *.o runStopDiag0L Root/*.o util/*.o


