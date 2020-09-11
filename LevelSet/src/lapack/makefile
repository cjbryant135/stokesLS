
# Root directory
ROOTDIR = $(PWD)/../../

# Compiler names and flags
CPP = ${CXX}

INCLUDE = -I$(ROOTDIR) $(XINCLUDE) 

depend:
	for i in *.cpp ; do $(CPP) -E -MMD $(INCLUDE) $$i > /dev/null ;  done 


