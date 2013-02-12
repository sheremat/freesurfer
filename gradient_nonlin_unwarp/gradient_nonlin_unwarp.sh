#!/bin/bash -p

# wrapper around MATLAB functions for unwarping images with gradient
# nonlinearity geometric distortion

# jonathan polimeni <jonp@nmr.mgh.harvard.edu>
# Tuesday, January 25, 2011 17:23:40 -0500
# Monday, September 24, 2012  0:56:03 -0400

# future options:
# - write out image data file with magnitudes of displacements
# - precalulated table of evaluated spharm / direct spharm evaluation / full gradient table from model

# todo:
# - check if input file exists
# - confirm that gradient file name is supported
# - slurp in previously-computed displacement map for re-running correction on several identical volumes

# to apply shiftmap:
#
#  applywarp \
#      --ref=${instem}.nii \
#      --in=${instem}.nii \
#      --warp=${outstem}__warp_ABS.nii.gz \
#      --abs \
#      --interp sinc \
#      --out=${outstem}__apply.nii.gz \
#      --datatype=float

# $Id: gradient_nonlin_unwarp.sh,v 1.6 2012/12/22 00:20:16 jonp Exp $
#**************************************************************************#

# set directory containing "mris_gradient_nonlin*.m" files
if [ -z ${MRIS} ]; then
    GRADNONLINUNWARP=${FREESURFER_HOME}/gradient_nonlin_unwarp
else
    GRADNONLINUNWARP=${MRIS}
fi
SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


#==========================================================================#

function askusage() {
    echo `basename ${0}` '<infile>  <outfile>  <coeffs>  [options]'
}

function askversion() {
    echo ' version 0.8   [$Revision: 1.6 $, $Date: 2012/12/22 00:20:16 $]'
}

function askhelp() {
    askusage `basename ${0}`
    echo -e ''
    echo -e ' ARGUMENTS'
    echo -e ''
    echo -e '   infile'
    echo -e '   outfile'
    echo -e ''
    echo -e '   coeffs           can be the path to a Siemens coeff.grad file or '
    echo -e '                    the name of a gradient coil on file (e.g., AC84, Bay4, connectome)'
    echo -e ''
    echo -e ' OPTIONS'
    echo -e ''
    echo -e '  --vol'
    echo -e '  --surf'
    echo -e ''
    echo -e '  --reg <regfile.dat>'
    echo -e '  --method <direct|lookup>'
    echo -e '  --polarity <UNDIS|DIS>'
    echo -e ''
    echo -e '  --biascor'
    echo -e '  --nobiascor'
    echo -e ''
    echo -e '  --interp <cubic|nearest|linear|spline>'
    echo -e ''
    echo -e '  --shiftmap'
    echo -e '  --noshiftmap'
    echo -e ''
    echo -e '  --savejac'
    echo -e ''
    echo -e '  --DEBUG #'
    echo -e ''
    echo -e ''
    echo -e ' jonathan polimeni <jonp@nmr.mgh.harvard.edu>'
    askversion
}

if [ ! ${#@} -gt 0 ]; then
    askusage `basename ${0}`
    exit 1
fi

if [ ${#@} -eq 1 ]; then
    case $1 in
        -h | -? | -help | --help)
            askhelp
            ;;
        -v | -version | --version)
            askversion
            ;;
    esac
    exit 1
fi


#--------------------------------------------------------------------------#

scriptname=$0
infile=$1
outfile=$2
coefffile=$3

# check to make sure that none of the required arguments are mis-interpreted
# as optional flags
if [[ "${infile}" == -* ]]; then
    askusage
    exit 1
fi
if [[ "${outfile}" == -* ]]; then
    askusage
    exit 1
fi
if [[ "${coefffile}" == -* ]]; then
    askusage
    exit 1
fi

cmd1="$*"
cmd0="$@"

shift 3

regfile=
template=
FLAG__surface=0
method='direct'
polarity='UNDIS'
biascor=1
interp='cubic'
JacDet=0
shiftmap=0
overlay=0
debug=0

for OPTION in "$@"; do
    case $OPTION in
        --surf)
            FLAG__surface=1
            shift 1
            ;;
        --vol)
            FLAG__surface=0
            shift 1
            ;;
        --reg)
            regfile=$2
            shift 2
            ;;
        --template)
            template=$2
            shift 2
            ;;
        --method)
            method=$2
            shift 2
            ;;
        --polarity)
            polarity=$2
            shift 2
            ;;
        --biascor)
            biascor=1
            echo "bias correction ENABLED"
            shift 1
            ;;
        --nobiascor)
            biascor=0
            echo "bias correction disabled"
            shift 1
            ;;
        --interp)
            interp=$2
            echo interp: $interp
            shift 2
            ;;
        --nearest)
            interp="nearest"
            echo interp: $interp
            shift 1
            ;;
        --linear)
            interp="linear"
            echo interp: $interp
            shift 1
            ;;
        --spline)
            interp="spline"
            echo interp: $interp
            shift 1
            ;;
        --cubic)
            interp="cubic"
            echo interp: $interp
            shift 1
            ;;
        --shiftmap)
            shiftmap=1
            echo "shift map output enabled"
            shift 1
            ;;
        --noshiftmap)
            shiftmap=0
            shift 1
            ;;
        --savejac)
            JacDet=1
            echo "saving jacobian determinant"
            shift 1
            ;;
        --debug)
            debug=1
            echo debug: $debug
            shift 1
            ;;
        [[:graph:]])
            echo "ERROR: option ""$1"" unrecognized"
            exit 1
            ;;
    esac
done


#--------------------------------------------------------------------------#

if (( $FLAG__surface == 1 )); then
    echo "surface warping selected..."
    cmd='mris_gradient_nonlin__unwarp_surface__batchmode'
else
    echo "volume warping selected..."
    cmd='mris_gradient_nonlin__unwarp_volume__batchmode'
fi


matlabcmd=`which matlab`

if [ ! -n "$matlabcmd" ]; then
    echo "ERROR: cannot find matlab in path."
    exit 1
fi


#--------------------------------------------------------------------------#

outstem=`basename ${outfile} .mgz | sed - -e s/\.nii\.gz$// | sed - -e s/\.nii$// | sed - -e s/\.mgh$//`

LOGFILE="`dirname ${outfile}`/${outstem}__`basename ${scriptname} .sh`.log"
echo LOGGING to file \"`basename $LOGFILE`\"

echo "$scriptname $cmd0" > ${LOGFILE}


if (( $FLAG__surface == 1 )); then

    echo -e "                ${cmd}('${infile}', '${outfile}', '${coefffile}', '${polarity}', '${method}', '${regfile}', '${template}', ${overlay}, ${debug})\n\n" | tee -a ${LOGFILE}

    ${matlabcmd} -nosplash -nodesktop \
        -r "addpath ${GRADNONLINUNWARP}; ${cmd}('${infile}', '${outfile}', '${coefffile}', '${polarity}', '${method}', '${regfile}', '${template}', ${overlay}, ${debug})" \
        -logfile ${LOGFILE}
else

    echo -e "                ${cmd}('${infile}', '${outfile}', '${coefffile}', '${polarity}', '${method}', '${biascor}', '${interp}', '${JacDet}', '${shiftmap}')\n\n" | tee -a ${LOGFILE}

    ${matlabcmd} -nosplash -nodesktop \
        -r "addpath ${GRADNONLINUNWARP}; ${cmd}('${infile}', '${outfile}', '${coefffile}', '${polarity}', '${method}', '${biascor}', '${interp}', '${JacDet}', '${shiftmap}')" \
        -logfile ${LOGFILE}
fi


exit 0


#**************************************************************************#
# $Source: /usr/fscvsroot/dev/gradient_nonlin_unwarp/gradient_nonlin_unwarp.sh,v $
# Local Variables:
# mode: sh
# fill-column: 76
# comment-column: 0
# End:
