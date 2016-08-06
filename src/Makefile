#
# Makefile for Saturation project.
#
TESTS 	= testtranslation.exe		\
			saturationdatasizetest.exe


TARGETS = saturation.exe 		\
		canonoutput.exe	


OBJECTS = 	Augmentation.o				\
			connectivity.o				\
			graphcompactor.o			\
			SaturationAugmenter.o		\
			SaturationData.o			\
			SaturationDataSymmetry.o	\
			SaturationGraph.o 			\
			SaturationManager.o 		\
			SaturationSymmetry.o		\
			Set.o						\
			translation.o				\
			TreeSet.o	


CC 	= gcc
CXX = g++

OPT = -O4
DEBUG = 
WARNINGS = -Wall -Werror

CFLAGS 		= $(OPT) $(DEBUG) $(WARNINGS)
CXXFLAGS 	= $(OPT) $(DEBUG) $(WARNINGS)
LFLAGS 		= $(OPT) $(DEBUG) $(WARNINGS)

TREESEARCHDIR = ../../TreeSearch
NAUTYDIR = ../../nauty

INCLUDES =  -I$(TREESEARCHDIR) \
			-I$(NAUTYDIR) 

TREESEARCHOBJS	= $(TREESEARCHDIR)/SearchManager.o 

NAUTYOBJS     	= $(NAUTYDIR)/nauty.o 			\
				  $(NAUTYDIR)/nausparse.o		\
				  $(NAUTYDIR)/gtools.o			\
				  $(NAUTYDIR)/nautil.o			\
				  $(NAUTYDIR)/naugraph.o		\
				  $(NAUTYDIR)/schreier.o		\
				  $(NAUTYDIR)/naurng.o

LIBOBJS			= $(TREESEARCHOBJS) $(NAUTYOBJS)

.SUFFIXES: .c .cpp .o .obj .exe 

all: $(OBJECTS) $(TESTS) $(TARGETS)

.c.o: $<
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -c $< -o $@
        
.cpp.o: $<
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LIBS) -c $< -o $@
        
.cpp.exe: $< $(OBJECTS)  
	$(CXX) $(LFLAGS)			\
	$(INCLUDES)	$(DEBUG)			\
	$(LIBOBJS) $(LIBS)				\
	`cat $@.objs`           		\
	$< -o $@

.c.exe: $< $(COBJECTS)
	$(CC) 	$(LFLAGS)			    \
	$(INCLUDES)				\
	$(NAUTYOBJS)  $(COBJECTS) $(LIBS)		\
	$< -o $@

clean:
	-@rm $(OBJECTS) $(TARGETS) $(TESTS)

cleanexe:
	-@rm $(TARGETS)

cleantest:
	-@rm $(TESTS)

clexe:
	-@rm $(TARGETS)
