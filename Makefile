TARGET = TMBResidualMain 
CC = g++
ROOTCINT = $(ROOTSYS)/bin/rootcint
DICTNAME = TMBResidualClass_dict
SRCS = $(addsuffix .C, $(TARGET))
DIR = .
SOURCES = $(DIR)/TMBResidualClass.cc
HEADERS = $(DIR)/TMBResidualClass.hh
#DEPS = $(DIR)/reconmcvar.hh $(DIR)/vetoreadcfg.hh
OBJS = $(addsuffix .o, $(notdir $(basename $(SRCS))))
SOBJS = $(addsuffix .o, $(notdir $(basename $(SOURCES))))
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags) -std=c++0x #-Wall -fPIC -g
ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs) -lProof -lProofPlayer #-lRooFit -lRooFitCore
ROOTLDFLAGS = $(shell $(ROOTSYS)/bin/root-config --ldflags)

all: $(DICTNAME).C
	$(CC) $(ROOTCFLAGS) $(ROOTLIBS) $(SRCS) $(SOURCES) $^ -o $(TARGET) 

$(DICTNAME).C: $(HEADERS) $(DEPS)
	$(ROOTCINT) -f $@ -c  $^ 

.PHONY: clean

clean:
	rm -rf $(DIR)/*_dict.h $(DIR)/*_dict.C