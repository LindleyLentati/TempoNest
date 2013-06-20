
ACLOCAL_AMFLAGS = -I autoconf

# Source files used by plugins that are not part of temponest
xtra_dist_sources=triplot.py

# Specify to include the T2runtime in distributions
EXTRA_DIST = $(xtra_dist_sources)

# This clears out any CVS files that might get in our distribution
#dist-hook:
#	rm -rf `find $(distdir)/T2runtime -name CVS`
#	rm -rf $(distdir)/T2runtime/plugins/*
#	rm -rf $(distdir)/unsupported_plugins/built_plugins/*
#	rm -rf `find $(distdir)/unsupported_plugins/ -name CVS`
#	rm -rf `find $(distdir)/unsupported_plugins/ -name '*.o'`



plugindir=@TEMPO2_PLUG_INSTALL@

if HAVE_MULTINEST

plugin_LTLIBRARIES=temponest_@T2ARCH@_plug.la

if HAVE_CULA

temponest_@T2ARCH@_plug_la_SOURCES= tempo2pred_int.h T2accel.h tempo2pred.h tempo2Util.h constraints.h TKmatrix.h TKfit.h T2toolkit.h tempo2.h config.h configfile.h configfile.cpp TempoNestGPUFuncs.cu dpotrs.cpp dgesvd.cpp dgemm.cpp dgemv.cpp dpotri.cpp dpotrf.cpp TempoNestParams.cpp MultiNestParams.cpp TempoNestTextOutput.cpp TempoNestUtilities.cpp TempoNestUpdateLinear.cpp TempoNestFindMaxGPU.cpp TempoNestSim.cpp TempoNestLikeFuncs.cpp TempoNestGPULikeFuncs.cpp TempoNestPlugin.cpp
temponest_@T2ARCH@_plug_la_LIBADD=-ltempo2 -lsofa @MULTINEST_LIBS@ @CULA_LIBS@ @GSL_LIBS@ @BLAS_LIBS@ @LIBS_LAPACK@ @TEMPO2_LIBS@
temponest_@T2ARCH@_plug_la_INCLUDES=@MULTINEST_CFLAGS@ @CULA_CFLAGS@ @GSL_CFLAGS@ @TEMPO2_CFLAGS@

else

temponest_@T2ARCH@_plug_la_SOURCES= tempo2pred_int.h T2accel.h tempo2pred.h tempo2Util.h constraints.h TKmatrix.h TKfit.h T2toolkit.h tempo2.h config.h configfile.h configfile.cpp dgemm.cpp dgesvd.h dpotrs.cpp TempoNestFindMax.cpp TempoNest.h dgemm.h dpotrf.cpp dpotrs.h TempoNestLikeFuncs.cpp TempoNestSim.cpp dgemv.cpp dpotrf.h multinest.h TempoNestTextOutput.cpp dgemv.h dpotri.cpp MultiNestParams.cpp  TempoNestUpdateLinear.cpp dgesvd.cpp  dpotri.h TempoNestPlugin.cpp TempoNestParams.cpp TempoNestUtilities.cpp
temponest_@T2ARCH@_plug_la_LIBADD=-ltempo2 -lsofa @MULTINEST_LIBS@ @GSL_LIBS@ @BLAS_LIBS@ @LIBS_LAPACK@ @TEMPO2_LIBS@
temponest_@T2ARCH@_plug_la_INCLUDES=@MULTINEST_CFLAGS@ @GSL_CFLAGS@ @TEMPO2_CFLAGS@

endif


###
#
# This sets the includes and options that are needed to build nice modules with libtool
#
###
LDFLAGS+= -shared -module -avoid-version -shrext .t2
INCLUDES=-I@top_srcdir@/



temponest: temponest_@T2ARCH@_plug.la
temponest-install: temponest_@T2ARCH@_plug.la plugdir
	$(INSTALL) .libs/temponest_@T2ARCH@_plug.t2 $(plugindir)

plugdir:
	mkdir -p $(plugindir)


endif


#.PHONY: plugins-install plugins plugins-clean \
#   unsupported unsupported-install unsupported-clean \
#   complete complete-install complete-clean \
#   core core-install core-clean \
#   temponest temponest-install temponest-clean
#
#install-data-hook:
#	@echo ""
#	@echo ""
#	@echo "*** NOTICE ***"
#	@echo "Plugins are no longer installed by default"
#	@echo "Use: 'make plugins-install' to install them"
#	@echo ""
#	@echo ""
