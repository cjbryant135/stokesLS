#************************************************
#   makefile
#
#   $Header: makefile,v 1.6 99/01/04 13:22:31 chopp Exp $
#
#   $Log:	makefile,v $
# Revision 1.6  99/01/04  13:22:31  13:22:31  chopp (David Chopp)
# *** empty log message ***
# 
# Revision 1.3  97/12/04  10:17:07  10:17:07  chopp (David Chopp)
# *** empty log message ***
# 
# Revision 1.2  97/12/03  15:00:26  15:00:26  chopp (David Chopp)
# Added spotless target
# 
# Revision 1.1  97/12/01  15:56:43  15:56:43  chopp (David Chopp)
# Initial revision
# 
#************************************************/

# Tunable parameters

# Root directory
ROOTDIR = $(PWD)

# Name of the library
LIBNAME = level

# End of tunable parameters

# Location of the library source code
SRCDIR = $(ROOTDIR)/src
HEADDIR = $(ROOTDIR)/$(LIBNAME)

# Where the library will be stored
LIBDIR = $(ROOTDIR)/lib

# Full name of the library
LIBFULL = $(LIBDIR)/lib$(LIBNAME).a

# Where the docs are stored
DOXDIR = $(ROOTDIR)/docs
DOXTOP = $(DOXDIR)/index.html

# Where the test programs are stored
TESTDIR = $(ROOTDIR)/test

lib: nothing
	mkdir -p lib
	cd $(SRCDIR); make lib

all: lib docs

docs: nothing
	mkdir -p docs
	cd $(HEADDIR); make docs

edit: nothing
	emacs -title `basename ${PWD}` -cr green $(HEADDIR)/*.h $(HEADDIR)/*.cc $(SRCDIR)/*.cpp $(EXAMDIR)/*.cpp

smalledit: nothing
	emacs -title `basename ${PWD}` -cr green -font 6x10 $(HEADDIR)/*.h $(HEADDIR)/*.cc $(SRCDIR)/*.cpp $(EXAMDIR)/*.cpp

clean: nothing
	cd $(SRCDIR); make clean
	cd $(HEADDIR); make clean
	rm -f *~

# Make sure LABEL is set
freeze: nothing
	@-if [[ "${LABEL}" != "" ]] ; then                             \
	  cd $(SRCDIR); make freeze LABEL=${LABEL} ; cd $(ROOTDIR)  ; \
	  cd $(HEADDIR); make freeze LABEL=${LABEL} ; cd $(ROOTDIR) ; \
	  cd $(EXAMDIR); make freeze LABEL=${LABEL} ; cd $(ROOTDIR) ; \
	  cd $(DOXDIR); make freeze LABEL=${LABEL} ; cd $(ROOTDIR)  ; \
	  cd $(TESTDIR); make freeze LABEL=${LABEL} ; cd $(ROOTDIR) ; \
	  rcs_freeze ${LABEL} makefile                              ; \
	else                                                          \
	  echo "Error: LABEL must be set for freeze"               ; \
	fi

spotless: clean
	rm -rf docs;
	rm -f lib/liblevel*;
	rmdir lib;

nothing:

