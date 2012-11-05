# TARGETS = stdvectorDict.cxx generator

# alldicts += stdvectorDict.cxx

# include Rules.mk

# generator:	main.cc TwoBodyDecayGen.cxx
# generator:	LDLIBS += $(ROOTLIBS)

ROOTINCLUDE	   := $(shell root-config --incdir)
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)
ROOFITLIBS := -lRooFitCore -lRooFit

.PHONY.:	clean

# Use -std=c++11 if using new initialiser
generator:	stdvectorDict.cxx main.cc TwoBodyDecayGen.cxx
	g++ -Wall -g $(ROOTCFLAGS) $(ROOTLIBS) $(ROOFITLIBS) $^ -o $@

stdvectorDict.cxx:	stdvectorLinkDef.h
	rootcint -f $@ -c -p $^
	grep include $^ > tmpdict && cat $@ >> tmpdict && mv tmpdict $@

test:	stdvectorDict.cxx test.cc
	g++ -Wall -g $(ROOTCFLAGS) $(ROOTLIBS) $^ -o $@

testpartial:	stdvectorDict.cxx testpartial.cc
	g++ -Wall -g $(ROOTCFLAGS) $(ROOTLIBS) $^ -o $@

clean:
	rm -rf generator stdvectorDict.{h,cxx} test testpartial

# libPhaseSpaceGen.so:
# rootcint -f ${TGT} -c -p ${CXXFLAGS} ${ROOTCINT_INCLUDES} ${SRC};
#         'grep include ${SRC} > ${TGT}-tmp || true; '
#         'cat ${TGT} >> ${TGT}-tmp; '
#         'mv ${TGT}-tmp ${TGT}
