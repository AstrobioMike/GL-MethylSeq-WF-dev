#!/usr/bin/env bash

# from: https://gist.github.com/toniher/65ccf76a8903bf432435d490b2025fab
# Addresses issue: https://github.com/nextflow-io/nextflow/issues/1210

# base usage: bash prepull_singularity.sh
    # this will look for an NXF_SINGULARITY_CACHEDIR environment variable and store them there,
    # if that doesn't exist, it will put them in a subdirectory called "singularity" in the current working directory,
    # and download the containers listed in config/software/docker-images.config by default (suitable for the workflow this is packaged with)

# if wanting to specify the storage location and configuration file
    # the first positional argument is the storage location
    # the second is the location of the configuration file, e.g.: bash prepull_singularity.sh my-containers/ different-configuration.config


if [ ! -z ${NXF_SINGULARITY_CACHEDIR} ]; then

    OUTDIR=${NXF_SINGULARITY_CACHEDIR}

else

    OUTDIR=${1:-./singularity}

fi

CONFILE=${2:-config/software/docker-images.config}

if [ ! -e ${CONFILE} ]; then
    echo "${CONFILE} does not exist"
    exit
fi

printf "\n    Configuarion file being used:\n        ${CONFILE}\n\n"
printf "    Location of stored images:\n        ${OUTDIR}\n\n"

TMPFILE=$(mktemp)

CURDIR=$(pwd)

mkdir -p $OUTDIR

cat ${CONFILE} | grep 'container' | perl -lane 'if ( $_=~/container\s*\=\s*\"(\S+)\"/ ) { $_=~/container\s*\=\s*\"(\S+)\"/; print $1 unless ( $1=~/^\s*$/ or $1=~/\.sif/ or $1=~/\.img/ ) ; }' > $TMPFILE

cd ${OUTDIR}

while IFS= read -r line; do
    name=$line
    name=${name/:/-}
    name=${name//\//-}
    echo $name
    singularity pull ${name}.img docker://$line
done < $TMPFILE

cd $CURDIR
