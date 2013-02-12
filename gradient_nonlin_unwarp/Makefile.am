#
# Makefile.am
#

SUBDIRS=gradient_coil_files

unwarpdir=$(prefix)/gradient_nonlin_unwarp

if NMR_INSTALL
unwarp_DATA=\
	mris_gradient_nonlin__jacobian_determinant.m \
	mris_gradient_nonlin__load_gradtable.m \
	mris_gradient_nonlin__read_siemens_coeff.m \
	mris_gradient_nonlin__setup.m \
	mris_gradient_nonlin__spharm_evaluate.m \
	mris_gradient_nonlin__unwarp.m \
	mris_gradient_nonlin__unwarp_surface__batchmode.m \
	mris_gradient_nonlin__unwarp_surface.m \
	mris_gradient_nonlin__unwarp_volume__batchmode.m \
	mris_gradient_nonlin__unwarp_volume.m \
	mris_read.m \
	mris_read_mgh.m \
	mris_read_mgh_vox2tkr.m \
	mris_read_nii.m \
	mris_read_surface.m \
	mris_read_vox2tkr.m \
	mris_save_surface.m \
	mris_transform_coordinates.m \
	mris_transform.m \
	mris_read_surface_fs.m \
	mris_read_surface_volgeom.m \
	mris_save_surface_fs.m \
	mris_show_matrix.m \
	timing.m \
	gradient_nonlin_unwarp.sh
else
unwarp_DATA=
endif

EXTRA_DIST=$(unwarp_DATA)

if NMR_INSTALL
install-data-hook:
	chmod a+x $(DESTDIR)$(unwarpdir)/gradient_nonlin_unwarp.sh
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=$(unwarp_DATA)
include $(top_srcdir)/Makefile.extra
