
# Root directory
ROOTDIR = $(PWD)/../../

# Compiler names and flags
CPP = ${CXX}
DEPEND = -DAVS_WRITE_UNIFORM -DPERIODIC_BOUNDARY -DNETCDF
DBGFLAGS = -g -DLEVEL_DEBUG -DCHECK_DIVIDE_BY_ZERO -DUSE_INLINES
OPTFLAGS = -O3
CPPFLAGS = $(DBGFLAGS) $(DEPEND) 
#CPPFLAGS = -m64 $(OPTFLAGS) $(DEPEND)

INCLUDE = -I$(ROOTDIR) $(XINCLUDE) -I$(ROOTDIR)/level -I$(ROOTDIR)/level/initfuncs

depend:
	for i in *.cpp ; do $(CPP) -E -MMD $(INCLUDE) $$i > /dev/null ;  done 


objects:
	for i in *.cpp ; do $(CPP) $(CPPFLAGS) $(INCLUDE) -c $$i; done

clean:
	@- rm -f *.o *.dat $(TARGETS) *.dox *~ \#*\#
