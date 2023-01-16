#!/bin/bash

VERSION=`cat VERSION.txt`

docker build -t brianjohnhaas/evidencemodeler:${VERSION} .
docker build -t brianjohnhaas/evidencemodeler:latest .

