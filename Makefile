TARGETS = stdvectorDict.cxx libDecayGen.so generator test testpartial

alldicts += stdvectorDict.cxx

LIBSRC = TwoBodyDecayGen.cxx $(alldicts)
BINSRC = generator.cc test.cc testpartial.cc

include mk/Rules.mk

# Common
$(TARGETS):	LDLIBS += -lstdc++ $(ROOTLIBS)

# Libraries
stdvectorDict.cxx:	stdvectorInclude.h stdvectorLinkDef.h

libDecayGen.so: $(patsubst %.cxx,%.os,$(filter-out %.cc,$(ccsrc)))


# Binaries
generator:	LDLIBS += -L./ -lDecayGen

test:		LDLIBS += -L./ -lDecayGen

testpartial:	LDLIBS += -L./ -lDecayGen


# Documentation
.PHONY:	docs gh-pages

docs:
	mkdir -p docs
	doxygen doxy.conf > /dev/null

gh-pages:	docs
	cd docs && git add html && \
	           git commit -a -m "Update HTML docs" && \
	           git push



##
## Old Makefile contents
##

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
