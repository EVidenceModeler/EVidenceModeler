#!/bin/bash

VERSION=`cat VERSION.txt`

singularity build EVidenceModeler.v${VERSION}.simg docker://brianjohnhaas/evidencemodeler:$VERSION

singularity exec -e EVidenceModeler.v${VERSION}.simg EVidenceModeler --version
