bin_PROGRAMS = pgrt

pgrt_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

pgrt_CPPFLAGS = -I$(top_srcdir)/src \
	-I$(top_srcdir)/vendor \
	-I$(top_srcdir)

pgrt_LDFLAGS = $(OPENMP_CXXFLAGS) \
	-pthread

pgrt_LDADD = -lz
	
pgrt_SOURCES = PanGenomeRefBuild.cpp \
	SingleGenomeUnique.hpp \
	Options.h \
	Util.h
	