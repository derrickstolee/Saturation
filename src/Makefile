#
# Base Makefile for EarSearch project.
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

<<<<<<< HEAD
OPT = -O4 -static
# -O7
=======
OPT = -O4 
>>>>>>> 54e39aa63d051762fc09f22b848b5047ef0e72be
DEBUG =  
# -g 
# -static 
WARNINGS = -Wall


CFLAGS 		= $(OPT) $(DEBUG) $(WARNINGS)
CXXFLAGS 	= $(OPT) $(DEBUG) $(WARNINGS)
LFLAGS 		= $(OPT) $(DEBUG) $(WARNINGS)


INCLUDES =  -I../TreeSearch \
			-I../nauty 


TREESEARCHOBJS	= ../TreeSearch/SearchManager.o 
				  
NAUTYOBJS     	= ../nauty/nauty.o 			\
				  ../nauty/nausparse.o		\
				  ../nauty/gtools.o			\
				  ../nauty/nautil.o			\
				  ../nauty/naugraph.o
					
LIBOBJS			= $(TREESEARCHOBJS) $(NAUTYOBJS)
			
			
.SUFFIXES: .c .cpp .o .obj .exe 

all: $(OBJECTS) $(TESTS) $(TARGETS)


# The default object compiler
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
	rm $(OBJECTS) $(TARGETS) $(TESTS)
	
cleanexe:
	rm $(TARGETS)

cleantest:
	rm $(TESTS)

clexe:
	rm $(TARGETS)
