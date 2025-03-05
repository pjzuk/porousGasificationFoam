#!/bin/bash

foamCleanPolyMesh
rm -r 0
cp -r save  0
cp initialData/radiationProperties_meantemp constant/radiationProperties

blockMesh
setFields
renumberMesh -overwrite
checkMesh
