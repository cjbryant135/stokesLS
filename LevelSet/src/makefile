#*************************************************
#   makefile
#
#   $Header: makefile,v 2.3 99/01/06 13:59:29 chopp Exp $
#
#   $Log:	makefile,v $
# Revision 2.3  99/01/06  13:59:29  13:59:29  chopp (David Chopp)
# *** none ***
# 
# Revision 1.8  98/05/18  14:06:30  14:06:30  chopp (David Chopp)
# *** empty log message ***
# 
# Revision 1.7  98/03/04  14:59:51  14:59:51  chopp (David Chopp)
# *** empty log message ***
# 
# Revision 1.7  98/03/02  12:58:18  12:58:18  chopp (David Chopp)
# *** empty log message ***
# 
# Revision 1.6  97/12/04  10:07:20  10:07:20  chopp (David Chopp)
# Added targets for boundary condition code and patched up some dependency
# information
# 
# Revision 1.5  97/12/03  15:17:00  15:17:00  chopp (David Chopp)
# *** empty log message ***
# 
# Revision 1.4  97/12/03  14:21:24  14:21:24  chopp (David Chopp)
# Turned the templatized code into non-templatized code
# 
# Revision 1.3  97/12/03  14:03:21  14:03:21  chopp (David Chopp)
# Added cakeflux1.cpp and cakeflux2.cpp to RCSLIST
# 
# Revision 1.2  97/12/01  15:37:28  15:37:28  chopp (David Chopp)
# Made the make freeze command ignore errors
# 
# Revision 1.1  97/12/01  15:31:11  15:31:11  chopp (David Chopp)
# Initial revision
# 
#************************************************/


# Tunable parameters

# Root directory
ROOTDIR = $(PWD)/../

# Compiler names and flags
CPP = ${CXX}

#DEPEND = -DAVS_WRITE_UNIFORM -DPERIODIC_BOUNDARY -DNETCDF
#DEPEND = -DAVS_WRITE_UNIFORM -DPERIODIC_BOUNDARY -DHDF5
DEPEND = -DAVS_WRITE_UNIFORM -DPERIODIC_BOUNDARY
DBGFLAGS = -g -O0 -DLEVEL_DEBUG -DCHECK_DIVIDE_BY_ZERO -DUSE_INLINES
OPTFLAGS = -O3
CPPFLAGS = -m64 $(DBGFLAGS) $(DEPEND)\
					 -D_BSD_SOURCE -D_POSIX_SOURCE -D_FILE_OFFSET_BITS=64\
					 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE

# Archiver name
ARCOMMAND = ar -r

# Name of the library
LIBNAME = level

# End of tunable parameters

# Location of the library source code
SRCDIR = $(ROOTDIR)/src
HEADDIR = $(ROOTDIR)/$(LIBNAME)

# Where the library will be stored
LIBDIR = $(ROOTDIR)/lib

# Full name of the library
LIBFULL = $(LIBDIR)/lib$(LIBNAME)

# Where the docs are stored
DOXDIR = $(ROOTDIR)/docs
DOXTOP = $(DOXDIR)/index.html

# Compiler flags
INCLUDE = -I$(ROOTDIR) $(XINCLUDE) -I$(ROOTDIR)/level

# Source code list

SOURCES = boundary2d.cpp    heapelt2d.cpp    um2advect.cpp    um2boundary.cpp \
          um2deriv.cpp      um2heap.cpp      um2init.cpp      um2io.cpp       \
          uniformmesh2d.cpp mesh2d.cpp       interface.cpp    initialfunc.cpp \
                                                                              \
          um3advect.cpp     um3deriv.cpp     um3heap.cpp      um3init.cpp     \
          uniformmesh3d.cpp um3io.cpp        interface3d.cpp  heapelt3d.cpp   \
          initialfunc3d.cpp boundary3d.cpp   um3boundary.cpp                  \
                                                                              \
          spglobals.cpp     spmatrix.cpp     spqmatrix.cpp    blockvector.cpp \
          spvector.cpp      blockqmatrix.cpp blockglobals.cpp                 \
                                                                              \
          plotwindow.cpp    plotwindow2d.cpp plotwindow3d.cpp                 \
          multiplotwindow.cpp miscio.cpp                                      \
                                                                              \
          bitarray.cpp      bitmatrix.cpp    bitvector.cpp                    \
          debug.cpp         inputparams.cpp  segment.cpp      bicubic.cpp     \
          tricubic.cpp      bicubicgrid.cpp  tricubicgrid.cpp                 \
                                                                              \
          tinystr.cpp       tinyxml.cpp      tinyxmlerror.cpp                 \
          tinyxmlparser.cpp units.cpp                                         \
                                                                              \
          genbicubic.cpp    polarmesh2d.cpp  pm2io.cpp        pm2advect.cpp   \
          pm2init.cpp       pm2heap.cpp      pm2deriv.cpp     pm2boundary.cpp 

#          lapack/dgesvd.cpp
          

# SOURCES2 = initfuncs/circle.cpp     initfuncs/constant.cpp  \
#            initfuncs/diamond.cpp    initfuncs/heavifunc.cpp \
#            initfuncs/ifuncunion.cpp initfuncs/line.cpp      \
#            initfuncs/linearfunc.cpp initfuncs/oval.cpp      \
#            initfuncs/sine.cpp       initfuncs/square.cpp    \
#            initfuncs/wavyfunc.cpp   initfuncs/cup.cpp       \
#            initfuncs/dbloval.cpp    initfuncs/dblrect.cpp   \
#            initfuncs/hardsmile.cpp  initfuncs/arc.cpp       \
#            initfuncs/msquares.cpp   initfuncs/mdisks.cpp    \
#            initfuncs/plus.cpp       initfuncs/plus2.cpp     \
#            initfuncs/star.cpp                               \
#                                                                     \
#            um2linear.cpp  um2periodic.cpp  um2xperiodic.cpp um2yperiodic.cpp        \
#            boundaries/mirrorbc.cpp  boundaries/quarterbc.cpp        \
#                                                                     \
#            initfuncs/mound.cpp      initfuncs/ovoid.cpp     \
#            initfuncs/wobble.cpp     initfuncs/mcubes.cpp    \
#            initfuncs/mballs.cpp     initfuncs/sheet.cpp     \
#            initfuncs/snake.cpp                    \
#                                                                     \
#            um3linear.cpp  um3xyperiodic.cpp um3periodic.cpp         

SOURCES2 = um2linear.cpp  um2periodic.cpp  um2xperiodic.cpp um2yperiodic.cpp        \
           um3linear.cpp  um3xyperiodic.cpp um3periodic.cpp pm2linear.cpp           

OBJECTS = ${SOURCES:.cpp=.o}
OBJECTS2 = ${SOURCES2:.cpp=.o}
DEPENDS = ${SOURCES:.cpp=.d} ${SOURCES2:.cpp=.d}
RCSLIST = ${SOURCES} ${SOURCES2}

.SUFFIXES: .cpp .o .d

lib: $(LIBFULL).a $(LIBFULL)init.a
debuglib: $(LIBFULL)g.a
parlib: $(LIBFULL)p.a
initlib: $(LIBFULL)init.a

$(LIBFULL).a: $(OBJECTS)
	-if [ ! -d $(LIBDIR) ] ; then mkdir $(LIBDIR) ; fi
	$(ARCOMMAND) $(LIBFULL).a $?
	if [ -f /usr/bin/ranlib -o -f /bin/ranlib ] ; \
            then ranlib $(LIBFULL).a ; fi

$(LIBFULL)g.a: $(OBJECTS)
	-if [ ! -d $(LIBDIR) ] ; then mkdir $(LIBDIR) ; fi
	$(ARCOMMAND) $(LIBFULL)g.a $?
	if [ -f /usr/bin/ranlib -o -f /bin/ranlib ] ; \
            then ranlib $(LIBFULL)g.a ; fi

$(LIBFULL)p.a: $(OBJECTS)
	-if [ ! -d $(LIBDIR) ] ; then mkdir $(LIBDIR) ; fi
	$(ARCOMMAND) $(LIBFULL)p.a $?
	if [ -f /usr/bin/ranlib -o -f /bin/ranlib ] ; \
            then ranlib $(LIBFULL)p.a ; fi

$(LIBFULL)init.a: $(OBJECTS2) libdirs
	-if [ ! -d $(LIBDIR) ] ; then mkdir $(LIBDIR) ; fi
	$(ARCOMMAND) $(LIBFULL)init.a $(OBJECTS2) `ls initfuncs/*.o` `ls boundaries/*.o`
	if [ -f /usr/bin/ranlib -o -f /bin/ranlib ] ; \
            then ranlib $(LIBFULL)init.a ; fi

libdirs: 
	cd initfuncs; make objects; cd ..
	cd boundaries; make objects; cd ..

depend: $(DEPENDS)
	-rm make.depend;
	cd initfuncs; make depend; cd ..
	cd boundaries; make depend; cd ..
	cd lapack; make depend; cd ..
	cat *.d lapack/*.d initfuncs/*.d boundaries/*.d >> make.depend 
	rm *.d lapack/*.d initfuncs/*.d boundaries/*.d 

# Source code rules with dependencies

bitarray.o: bitarray.cpp ../level/bitarray.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

bitmatrix.o: bitmatrix.cpp ../level/bitmatrix.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

bitvector.o: bitvector.cpp ../level/bitvector.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

boundary2d.o: boundary2d.cpp ../level/boundary2d.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

cake.o: cake.cpp ../level/cake.h ../level/plotwindow2d.h \
 ../level/plotwindow.h ../level/initialfunc.h ../level/defs.h \
 ../level/datavec.h ../level/datablock.h ../level/reference.h \
 ../level/datablock.cc ../level/datavec.cc ../level/mesh2d.h \
 ../level/compiler.h ../level/interface.h ../level/segment.h \
 ../level/utility.h ../level/numtrait.h \
 ../level/bitmatrix.h ../level/sparsemat.h \
 ../level/debug.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

cakeadvect.o: cakeadvect.cpp ../level/numtrait.h ../level/datavec.h \
 ../level/datablock.h ../level/reference.h ../level/defs.h \
 ../level/datablock.cc ../level/datavec.cc ../level/cake.h \
 ../level/plotwindow2d.h ../level/plotwindow.h ../level/initialfunc.h \
 ../level/mesh2d.h ../level/compiler.h ../level/interface.h \
 ../level/segment.h ../level/utility.h \
 ../level/bitmatrix.h ../level/sparsemat.h \
 ../level/cakeheap.h ../level/heapt.h ../level/heapt.cc \
 ../level/cakeheapelt.h ../level/heapelt.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

cakederiv.o: cakederiv.cpp ../level/cake.h ../level/plotwindow2d.h \
 ../level/plotwindow.h ../level/initialfunc.h ../level/defs.h \
 ../level/datavec.h ../level/datablock.h ../level/reference.h \
 ../level/datablock.cc ../level/datavec.cc ../level/mesh2d.h \
 ../level/compiler.h ../level/interface.h ../level/segment.h \
 ../level/utility.h ../level/numtrait.h \
 ../level/bitmatrix.h ../level/sparsemat.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

cakeflux.o: cakeflux.cpp ../level/numtrait.h ../level/datavec.h \
 ../level/datablock.h ../level/reference.h ../level/defs.h \
 ../level/datablock.cc ../level/datavec.cc ../level/cake.h \
 ../level/plotwindow2d.h ../level/plotwindow.h ../level/initialfunc.h \
 ../level/mesh2d.h ../level/compiler.h ../level/interface.h \
 ../level/segment.h ../level/utility.h \
 ../level/bitmatrix.h ../level/sparsemat.h \
 ../level/cakeheap.h ../level/heapt.h ../level/heapt.cc \
 ../level/cakeheapelt.h ../level/heapelt.h cakeflux1.cpp cakeflux2.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

cakeheap.o: cakeheap.cpp ../level/cakeheap.h ../level/cake.h \
 ../level/plotwindow2d.h ../level/plotwindow.h ../level/initialfunc.h \
 ../level/defs.h ../level/datavec.h ../level/datablock.h \
 ../level/reference.h ../level/datablock.cc ../level/datavec.cc \
 ../level/mesh2d.h ../level/compiler.h ../level/interface.h \
 ../level/segment.h ../level/utility.h ../level/numtrait.h \
 ../level/bitmatrix.h \
 ../level/sparsemat.h ../level/heapt.h ../level/heapt.cc \
 ../level/cakeheapelt.h ../level/heapelt.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

cakeinit.o: cakeinit.cpp ../level/defs.h ../level/cake.h \
 ../level/plotwindow2d.h ../level/plotwindow.h ../level/initialfunc.h \
 ../level/datavec.h ../level/datablock.h ../level/reference.h \
 ../level/datablock.cc ../level/datavec.cc ../level/mesh2d.h \
 ../level/compiler.h ../level/interface.h ../level/segment.h \
 ../level/utility.h ../level/numtrait.h \
 ../level/bitmatrix.h ../level/sparsemat.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

cakeio.o: cakeio.cpp ../level/defs.h ../level/cake.h \
 ../level/plotwindow2d.h ../level/plotwindow.h ../level/initialfunc.h \
 ../level/datavec.h ../level/datablock.h ../level/reference.h \
 ../level/datablock.cc ../level/datavec.cc ../level/mesh2d.h \
 ../level/compiler.h ../level/interface.h ../level/segment.h \
 ../level/utility.h ../level/numtrait.h \
 ../level/bitmatrix.h ../level/sparsemat.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

debug.o: debug.cpp ../level/debug.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

hardsmile.o: hardsmile.cpp ../level/hardsmile.h ../level/initialfunc.h \
 ../level/defs.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

star.o: star.cpp ../level/star.h ../level/initialfunc.h ../level/defs.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

mirroredcubic.o: mirroredcubic.cpp ../level/mirroredcubic.h ../level/initialfunc.h ../level/defs.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

heapelt.o: heapelt.cpp ../level/heapelt.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

initialfunc.o: initialfunc.cpp ../level/initialfunc.h ../level/defs.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

inputparams.o: inputparams.cpp ../level/inputparams.h \
 ../level/initialfunc.h ../level/defs.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

interface.o: interface.cpp ../level/interface.h ../level/datavec.h \
 ../level/datablock.h ../level/reference.h ../level/defs.h \
 ../level/datablock.cc ../level/datavec.cc ../level/plotwindow2d.h \
 ../level/plotwindow.h ../level/segment.h ../level/initialfunc.h \
 ../level/numtrait.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

mesh2d.o: mesh2d.cpp ../level/mesh2d.h ../level/compiler.h \
 ../level/initialfunc.h ../level/defs.h ../level/interface.h \
 ../level/datavec.h ../level/datablock.h ../level/reference.h \
 ../level/datablock.cc ../level/datavec.cc ../level/plotwindow2d.h \
 ../level/plotwindow.h ../level/segment.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

miscio.o: miscio.cpp ../level/miscio.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

multiplotwindow.o: multiplotwindow.cpp ../level/multiplotwindow.h \
 ../level/plotwindow.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

oval.o: oval.cpp ../level/oval.h ../level/initialfunc.h \
 ../level/defs.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

plotwindow.o: plotwindow.cpp ../level/plotwindow.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

plotwindow2d.o: plotwindow2d.cpp ../level/plotwindow2d.h \
 ../level/plotwindow.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

plotwindow3d.o: plotwindow3d.cpp ../level/plotwindow3d.h \
 ../level/plotwindow.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

segment.o: segment.cpp ../level/segment.h ../level/plotwindow2d.h \
 ../level/plotwindow.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

sparsemat.o: sparsemat.cpp ../level/sparsemat.h ../level/datavec.h \
 ../level/datablock.h ../level/reference.h ../level/defs.h \
 ../level/datablock.cc ../level/datavec.cc
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

tinyxml.o: tinyxml.cpp ../level/tinyxml.h ../level/tinystr.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

tinyxmlerror.o: tinyxmlerror.cpp ../level/tinystr.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

tinyxmlparser.o: tinyxmlparser.cpp ../level/tinystr.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

um2advect.o: um2advect.cpp ../level/numtrait.h ../level/datavec.h \
 ../level/datablock.h ../level/reference.h ../level/defs.h \
 ../level/datablock.cc ../level/datavec.cc ../level/uniformmesh2d.h \
 ../level/plotwindow2d.h ../level/plotwindow.h ../level/initialfunc.h \
 ../level/mesh2d.h ../level/compiler.h ../level/interface.h \
 ../level/segment.h ../level/utility.h \
 ../level/um2heap.h ../level/heapt.h \
 ../level/heapt.cc ../level/um2heapelt.h ../level/heapelt.h \
 ../level/debug.h ../level/um2boundary.h ../level/boundary2d.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

um2deriv.o: um2deriv.cpp ../level/uniformmesh2d.h \
 ../level/plotwindow2d.h ../level/plotwindow.h ../level/initialfunc.h \
 ../level/defs.h ../level/datavec.h ../level/datablock.h \
 ../level/reference.h ../level/datablock.cc ../level/datavec.cc \
 ../level/mesh2d.h ../level/compiler.h ../level/interface.h \
 ../level/segment.h ../level/utility.h ../level/numtrait.h \
 ../level/um2boundary.h \
 ../level/boundary2d.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

um2heap.o: um2heap.cpp ../level/um2heap.h ../level/uniformmesh2d.h \
 ../level/plotwindow2d.h ../level/plotwindow.h ../level/initialfunc.h \
 ../level/defs.h ../level/datavec.h ../level/datablock.h \
 ../level/reference.h ../level/datablock.cc ../level/datavec.cc \
 ../level/mesh2d.h ../level/compiler.h ../level/interface.h \
 ../level/segment.h ../level/utility.h ../level/numtrait.h \
 ../level/heapt.h \
 ../level/heapt.cc ../level/um2heapelt.h ../level/heapelt.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

um2init.o: um2init.cpp ../level/defs.h ../level/uniformmesh2d.h \
 ../level/plotwindow2d.h ../level/plotwindow.h ../level/initialfunc.h \
 ../level/datavec.h ../level/datablock.h ../level/reference.h \
 ../level/datablock.cc ../level/datavec.cc ../level/mesh2d.h \
 ../level/compiler.h ../level/interface.h ../level/segment.h \
 ../level/utility.h ../level/numtrait.h \
 ../level/um2boundary.h ../level/boundary2d.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

um2io.o: um2io.cpp ../level/defs.h ../level/uniformmesh2d.h \
 ../level/plotwindow2d.h ../level/plotwindow.h ../level/initialfunc.h \
 ../level/datavec.h ../level/datablock.h ../level/reference.h \
 ../level/datablock.cc ../level/datavec.cc ../level/mesh2d.h \
 ../level/compiler.h ../level/interface.h ../level/segment.h \
 ../level/utility.h ../level/numtrait.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

um2linear.o: um2linear.cpp ../level/um2linear.h ../level/um2boundary.h \
 ../level/boundary2d.h ../level/uniformmesh2d.h \
 ../level/plotwindow2d.h ../level/plotwindow.h ../level/initialfunc.h \
 ../level/defs.h ../level/datavec.h ../level/datablock.h \
 ../level/reference.h ../level/datablock.cc ../level/datavec.cc \
 ../level/mesh2d.h ../level/compiler.h ../level/interface.h \
 ../level/segment.h ../level/utility.h ../level/numtrait.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

um2boundary.o: um2boundary.cpp ../level/um2boundary.h \
 ../level/boundary2d.h ../level/uniformmesh2d.h \
 ../level/plotwindow2d.h ../level/plotwindow.h ../level/initialfunc.h \
 ../level/defs.h ../level/datavec.h ../level/datablock.h \
 ../level/reference.h ../level/datablock.cc ../level/datavec.cc \
 ../level/mesh2d.h ../level/compiler.h ../level/interface.h \
 ../level/segment.h ../level/utility.h ../level/numtrait.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

um2periodic.o: um2periodic.cpp ../level/um2periodic.h \
 ../level/um2boundary.h ../level/boundary2d.h ../level/uniformmesh2d.h \
 ../level/plotwindow2d.h ../level/plotwindow.h ../level/initialfunc.h \
 ../level/defs.h ../level/datavec.h ../level/datablock.h \
 ../level/reference.h ../level/datablock.cc ../level/datavec.cc \
 ../level/mesh2d.h ../level/compiler.h ../level/interface.h \
 ../level/segment.h ../level/utility.h ../level/numtrait.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

um2xperiodic.o: um2xperiodic.cpp ../level/um2xperiodic.h \
 ../level/um2boundary.h ../level/boundary2d.h ../level/uniformmesh2d.h \
 ../level/plotwindow2d.h ../level/plotwindow.h ../level/initialfunc.h \
 ../level/defs.h ../level/datavec.h ../level/datablock.h \
 ../level/reference.h ../level/datablock.cc ../level/datavec.cc \
 ../level/mesh2d.h ../level/compiler.h ../level/interface.h \
 ../level/segment.h ../level/utility.h ../level/numtrait.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

um2yperiodic.o: um2yperiodic.cpp ../level/um2yperiodic.h \
 ../level/um2boundary.h ../level/boundary2d.h ../level/uniformmesh2d.h \
 ../level/plotwindow2d.h ../level/plotwindow.h ../level/initialfunc.h \
 ../level/defs.h ../level/datavec.h ../level/datablock.h \
 ../level/reference.h ../level/datablock.cc ../level/datavec.cc \
 ../level/mesh2d.h ../level/compiler.h ../level/interface.h \
 ../level/segment.h ../level/utility.h ../level/numtrait.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

uniformmesh2d.o: uniformmesh2d.cpp ../level/uniformmesh2d.h \
 ../level/plotwindow2d.h ../level/plotwindow.h ../level/initialfunc.h \
 ../level/defs.h ../level/datavec.h ../level/datablock.h \
 ../level/reference.h ../level/datablock.cc ../level/datavec.cc \
 ../level/mesh2d.h ../level/compiler.h ../level/interface.h \
 ../level/segment.h ../level/utility.h ../level/numtrait.h \
 ../level/debug.h \
 ../level/um2boundary.h ../level/boundary2d.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

units.o: units.cpp ../level/units.h
	$(CPP) $(CPPFLAGS) -D_PATH_UNITSLIB='"$(ROOTDIR)/lib/units.dat"' $(INCLUDE) -c $*.cpp
	#$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

genbicubic.o: genbicubic.cpp \
 ../level/genbicubic.h  ../level/lapack.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

polarmesh2d.o: polarmesh2d.cpp  ../level/polarmesh2d.h \
 ../level/plotwindow2d.h  ../level/plotwindow.h  ../level/plotwindow3d.h \
 ../level/initialfunc.h  ../level/defs.h  ../level/inputparams.h \
 ../level/tinyxml.h  ../level/genbicubic.h  ../level/datavec.h \
 ../level/datablock.h  ../level/reference.h  ../level/datablock.cc \
 ../level/datavec.cc  ../level/utility.h  ../level/interface.h \
 ../level/segment.h  ../level/numtrait.h  ../level/debug.h \
 ../level/numtrait.h  ../level/debug.h  ../level/pm2boundary.h \
 ../level/boundary2d.h  ../level/polarmesh2d.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

pm2io.o: pm2io.cpp ../level/defs.h \
 ../level/polarmesh2d.h  ../level/plotwindow2d.h  ../level/plotwindow.h \
 ../level/plotwindow3d.h  ../level/initialfunc.h  ../level/defs.h \
 ../level/inputparams.h  ../level/tinyxml.h  ../level/genbicubic.h \
 ../level/datavec.h  ../level/datablock.h  ../level/reference.h \
 ../level/datablock.cc  ../level/datavec.cc  ../level/utility.h \
 ../level/interface.h  ../level/segment.h  ../level/numtrait.h \
 ../level/debug.h  ../level/pm2boundary.h  ../level/boundary2d.h \
 ../level/polarmesh2d.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

pm2advect.o:  pm2advect.cpp \
 ../level/numtrait.h  ../level/datavec.h  ../level/datablock.h \
 ../level/reference.h  ../level/defs.h  ../level/datablock.cc \
 ../level/datavec.cc  ../level/datavec.h  ../level/polarmesh2d.h \
 ../level/plotwindow2d.h  ../level/plotwindow.h  ../level/plotwindow3d.h \
 ../level/initialfunc.h  ../level/inputparams.h  ../level/tinyxml.h \
 ../level/genbicubic.h  ../level/utility.h  ../level/interface.h \
 ../level/segment.h  ../level/numtrait.h  ../level/debug.h \
 ../level/pm2heap.h  ../level/polarmesh2d.h  ../level/heapt.h \
 ../level/heapt.cc  ../level/pm2heapelt.h  ../level/heapelt2d.h \
 ../level/pm2boundary.h  ../level/boundary2d.h  ../level/utility.h \
 ../level/genbicubic.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

pm2init.o: pm2init.cpp \
 ../level/defs.h  ../level/polarmesh2d.h  ../level/plotwindow2d.h \
 ../level/plotwindow.h  ../level/plotwindow3d.h  ../level/initialfunc.h \
 ../level/defs.h  ../level/inputparams.h  ../level/tinyxml.h \
 ../level/genbicubic.h  ../level/datavec.h  ../level/datablock.h \
 ../level/reference.h  ../level/datablock.cc  ../level/datavec.cc \
 ../level/utility.h  ../level/interface.h  ../level/segment.h \
 ../level/numtrait.h  ../level/debug.h  ../level/pm2boundary.h \
 ../level/boundary2d.h  ../level/polarmesh2d.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

pm2heap.o: pm2heap.cpp \
 ../level/pm2heap.h  ../level/polarmesh2d.h  ../level/plotwindow2d.h \
 ../level/plotwindow.h  ../level/plotwindow3d.h  ../level/initialfunc.h \
 ../level/defs.h  ../level/inputparams.h  ../level/tinyxml.h \
 ../level/genbicubic.h  ../level/datavec.h  ../level/datablock.h \
 ../level/reference.h  ../level/datablock.cc  ../level/datavec.cc \
 ../level/utility.h  ../level/interface.h  ../level/segment.h \
 ../level/numtrait.h  ../level/debug.h  ../level/heapt.h \
 ../level/heapt.cc  ../level/pm2heapelt.h  ../level/heapelt2d.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

pm2deriv.o: pm2deriv.cpp \
 ../level/polarmesh2d.h  ../level/plotwindow2d.h  ../level/plotwindow.h \
 ../level/plotwindow3d.h  ../level/initialfunc.h  ../level/defs.h \
 ../level/inputparams.h  ../level/tinyxml.h  ../level/genbicubic.h \
 ../level/datavec.h  ../level/datablock.h  ../level/reference.h \
 ../level/datablock.cc  ../level/datavec.cc  ../level/utility.h \
 ../level/interface.h  ../level/segment.h  ../level/numtrait.h \
 ../level/debug.h  ../level/pm2boundary.h  ../level/boundary2d.h \
 ../level/polarmesh2d.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

pm2boundary.o: pm2boundary.cpp \
 ../level/pm2boundary.h  ../level/boundary2d.h  ../level/polarmesh2d.h \
 ../level/plotwindow2d.h  ../level/plotwindow.h  ../level/plotwindow3d.h \
 ../level/initialfunc.h  ../level/defs.h  ../level/inputparams.h \
 ../level/tinyxml.h  ../level/genbicubic.h  ../level/datavec.h \
 ../level/datablock.h  ../level/reference.h  ../level/datablock.cc \
 ../level/datavec.cc  ../level/utility.h  ../level/interface.h \
 ../level/segment.h  ../level/numtrait.h  ../level/debug.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $*.cpp

lapack/dgesvd.o: lapack/dgesvd.cpp /usr/local/rogue/rw/lapkdefs.h \
  /usr/local/rogue/rw/mathdefs.h /usr/local/rogue/rw/defs.h \
  /usr/local/rogue/rw/generic.h /usr/local/rogue/rw/mathwind.h \
  /usr/local/rogue/rw/lapkwind.h /usr/local/rogue/rw/bla.h \
  /usr/local/rogue/rw/dcomplex.h /usr/local/rogue/rw/lapack.h \
  /usr/local/rogue/rw/fortran.h
	$(CPP) $(CPPFLAGS) $(INCLUDE) -o lapack/dgesvd.o -c lapack/dgesvd.cpp

unittest: unittest.o
	$(CPP) -g $(CPPFLAGS) $(INCLUDE) -L../lib -o unittest unittest.o -llevel

# General rules without header file dependencies

.c.o:
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $< 

.cpp.o:
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $< 

.cpp.d:
	$(CPP) -E -MMD $(INCLUDE) $< > /dev/null

clean:
	@- rm -f *.o *.dat $(TARGETS) *.dox *~ \#*\#
	cd initfuncs; make clean; cd ..
	cd boundaries; make clean; cd ..

freeze:
	@-if [[ "${LABEL}" != "" ]] ; then           \
	  rcs_freeze ${LABEL} ${RCSLIST}           ; \
	  rcs_freeze ${LABEL} makefile             ; \
	else                                         \
	  print "Error: LABEL must be set"         ; \
	fi


nothing:
