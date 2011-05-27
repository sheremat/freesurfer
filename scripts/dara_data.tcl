##
## dara_data.tcl
##
## CVS Revision Info:
##    $Author: nicks $
##    $Date: 2011/03/02 00:04:35 $
##    $Revision: 1.3 $
##
## Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
##
## Terms and conditions for use, reproduction, distribution and contribution
## are found in the 'FreeSurfer Software License Agreement' contained
## in the file 'LICENSE' found in the FreeSurfer distribution, and here:
##
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
##
## Reporting: freesurfer@nmr.mgh.harvard.edu
##

LoadFunctionalOverlay /space/beast/1/users/dara/normalsnew/bold/edpanalmc/tal-ffx/allvfix sig /space/beast/1/users/dara/normalsnew/bold/edpanalmc/tal-ffx/register.data
Overlay_SetThreshold 2 5 1
SetCursor 0 96 67 128
SelectValuesByFuncVoxel 1
RedrawScreen
