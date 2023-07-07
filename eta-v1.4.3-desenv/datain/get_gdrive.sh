#!/bin/bash
fileid=$1
filename=$2
curl -Lb /tmp/gcokie "https://drive.google.com/uc?export=download&confirm=Uq6r&id=${fileid}" -o "${filename}"
exit 0
