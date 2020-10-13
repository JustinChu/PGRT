bin_PROGRAMS = pgrt-build

pgrt_build_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

pgrt_build_CPPFLAGS = -I$(top_srcdir)/src \
	-I$(top_srcdir)/vendor \
	-I$(top_srcdir)
		
pgrt_build_LDFLAGS = $(OPENMP_CXXFLAGS)

pgrt_build_SOURCES = PanGenomeRefBuild.cpp