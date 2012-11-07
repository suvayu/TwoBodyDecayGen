TARGETS = stdvectorDict.cxx libDecayGen.so generator test testpartial

alldicts += stdvectorDict.cxx

LIBSRC = TwoBodyDecayGen.cxx $(alldicts)
BINSRC = generator.cc test.cc testpartial.cc

include mk/Rules.mk

SRCS += $(LIBSRC)
SRCS += $(BINSRC)

stdvectorDict.cxx:	stdvectorLinkDef.h
	rootcint -f $@ -c -p $^
	grep include $^ > tmpdict && cat $@ >> tmpdict && mv tmpdict $@

libDecayGen.so: $(patsubst %.cxx,%.os,$(filter-out $(BINSRC),$(SRCS)))
libDecayGen.so: LDLIBS += -lstdc++ $(ROOTLIBS)


# Binaries
generator:	LDLIBS += -lstdc++ $(ROOTLIBS) -L./ -lDecayGen

test:		LDLIBS += -lstdc++ $(ROOTLIBS)

testpartial:	LDLIBS += -lstdc++ $(ROOTLIBS) -L./ -lDecayGen


# ROOTINCLUDE	   := $(shell root-config --incdir)
# ROOTCFLAGS := $(shell root-config --cflags)
# ROOTLIBS   := $(shell root-config --libs)
# ROOFITLIBS := -lRooFitCore -lRooFit

# .PHONY.:	clean

# # Use -std=c++11 if using new initialiser
# generator:	stdvectorDict.cxx main.cc TwoBodyDecayGen.cxx
# 	g++ -Wall -g $(ROOTCFLAGS) $(ROOTLIBS) $(ROOFITLIBS) $^ -o $@

# stdvectorDict.cxx:	stdvectorLinkDef.h
# 	rootcint -f $@ -c -p $^
# 	grep include $^ > tmpdict && cat $@ >> tmpdict && mv tmpdict $@

# test:	stdvectorDict.cxx test.cc
# 	g++ -Wall -g $(ROOTCFLAGS) $(ROOTLIBS) $^ -o $@

# testpartial:	stdvectorDict.cxx testpartial.cc
# 	g++ -Wall -g $(ROOTCFLAGS) $(ROOTLIBS) $^ -o $@

# clean:
# 	rm -rf generator stdvectorDict.{h,cxx} test testpartial
