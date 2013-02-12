#
# Makefile.am
#

foodir=$(prefix)/gradient_nonlin_unwarp/gradient_coil_files

if NMR_INSTALL
tables=\
	coeff_AC44.grad \
	coeff_AC84.grad \
	coeff_AC84.gwt \
	coeff_AC84.mat \
	coeff_AC88.grad \
	coeff_AS05.grad \
	coeff_AS092.grad \
	coeff_AS302.grad \
	coeff_GC99.grad \
	coeff_SC72.grad \
	coeff_Sonata.grad \
	mris_gradient_nonlin__list_coeff_files.m \
	mris_gradient_nonlin__pick_coeff_file.m
else
tables=
endif

foo_DATA=$(tables)

EXTRA_DIST=$(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=$(tables)
include $(top_srcdir)/Makefile.extra
