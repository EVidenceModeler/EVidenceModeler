#!/bin/bash

VERSION=`cat VERSION.txt`


simg_name="EVidenceModeler.v${VERSION}.simg"

if [ -e $simg_name ]; then
    rm -f $simg_name
fi

singularity build ${simg_name} docker://brianjohnhaas/evidencemodeler:$VERSION

singularity exec -e ${simg_name} EVidenceModeler --version

ln -sf ${simg_name} EVidenceModeler.latest.simg
