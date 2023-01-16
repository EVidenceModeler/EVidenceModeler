#!/bin/bash

VERSION=`cat VERSION.txt`

docker push brianjohnhaas/evidencemodeler:${VERSION} 
docker push brianjohnhaas/evidencemodeler:latest 
