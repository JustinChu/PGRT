bin_PROGRAMS = pgrt pgrtCount

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
	
pgrtCount_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

pgrtCount_CPPFLAGS = -I$(top_srcdir)/src \
	-I$(top_srcdir)/vendor \
	-I$(top_srcdir)

pgrtCount_LDFLAGS = $(OPENMP_CXXFLAGS) \
	-pthread

pgrtCount_LDADD = -lz
	
pgrtCount_SOURCES = CoverageEstmate.cpp \
	CountKmers.hpp \
	Options.h \
	Util.h
	