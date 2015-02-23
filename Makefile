include ../CGNS.options

#
# set architecture
#
ARCH = $(shell uname)

ifeq ($(ARCH),Linux)
OS = $(shell uname -i)
ifeq ($(OS),x86_64)
	ARCH = Linux-64
endif
endif

#COMP = g++
#LINK = g++
#CFLAGS = -c -O3
COMP = mpic++
LINK = mpic++
CFLAGS = -c -O3 -DPARALLEL
#CFLAGS = -c -g
DEFINES =

INCLUDE_PATH = -I./
LIBRARY_PATH = -L./

########################################################################
#  MAC OSX, using gcc compilers
########################################################################
ifeq ($(ARCH),Darwin)
	METIS_INCLUDE_PATH = /usr/local/metis/include
	METIS_LIB_PATH	   = /usr/local/metis/lib
	CGNS_INCLUDE_PATH  = /usr/local/include
	CGNS_LIB_PATH      = /usr/local/lib
	EXE_SUFFIX = MACOSX
endif

########################################################################
#  Linux-64, using intel compilers
########################################################################
ifeq ($(ARCH),Linux-64)
	METIS_INCLUDE_PATH = /usr/local/metis/include
	METIS_LIB_PATH	   = /usr/local/metis/lib
	CGNS_INCLUDE_PATH  = /simcenter/meshdev/cgns/cgns_3.0/release/3.0.5/src
	CGNS_LIB_PATH      = /simcenter/meshdev/cgns/cgns_3.0/release/3.0.5/src/LINUX64
	EXE_SUFFIX = LINUX64
endif

INCLUDE_PATH += -I$(METIS_INCLUDE_PATH)
LIBRARY_PATH += -L$(METIS_LIB_PATH)

INCLUDE_PATH += -I../UTIL -I../OCTREE_LIBRARY -I../SGIO
LIBRARY_PATH += -L../UTIL -L../OCTREE_LIBRARY -L../SGIO

LDFLAGS = -lmetis -loctree -lmesh_io -lUtility -lsgio

ifeq ($(HAS_CGNS),yes)
  DEFINES += -DHAVE_CGNS
  INCLUDE_PATH += -I$(CGNS_INCLUDE_PATH) -I../PCGNS
  LIBRARY_PATH += -L$(CGNS_LIB_PATH) -L../PCGNS
  LDFLAGS += -lpcgns -lcgns
endif


CNVC_SRCS = Conv.cpp \
	convert.cpp \
	create_hybrid_maps.cpp \
	split_tree_partition.cpp \
	merge.cpp \

CNVSRCS = $(CNVC_SRCS)
CNVOBJECTS = $(CNVC_SRCS:.cpp=.o)

Conv: mesh_io.a $(CNVOBJECTS)
	$(LINK) -o $@ $(CNVOBJECTS) $(LIBRARY_PATH) $(LDFLAGS)
	cp $@ $@.$(EXE_SUFFIX)
	cp $@ ../$@.$(EXE_SUFFIX)

D_SRCS=\
	mesh_io.cpp \
	StarCD.cpp \

DOBJECTS = $(D_SRCS:.cpp=.o)

mesh_io.a:$(DOBJECTS)
	ar rv lib$@ $(DOBJECTS)
	ranlib lib$@

install:
	cp Conv /simcenter/meshdev/Conv.$(EXE_SUFFIX)

clean:
	/bin/rm -f *.o
	/bin/rm -f libmesh_io.a

.cpp.o:
	$(COMP) $(INCLUDE_PATH) $(DEFINES) $(CFLAGS) $<

.c.o:
	$(COMP) $(INCLUDE_PATH) $(DEFINES) $(CFLAGS) $<

