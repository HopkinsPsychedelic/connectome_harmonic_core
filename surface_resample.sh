#!/bin/bash

subsurfdir=$1

refspheredir=$2






mris_convert ${subsurfdir}/lh.midthickness ${subsurfdir}/lh.midthickness.vtk

mris_convert ${subsurfdir}/rh.midthickness ${subsurfdir}/rh.midthickness.vtk

mris_convert ${subsurfdir}/lh.white ${subsurfdir}/lh.white.vtk

mris_convert ${subsurfdir}/rh.white ${subsurfdir}/rh.white.vtk

mris_convert ${subsurfdir}/lh.sphere.reg ${subsurfdir}/lh.sphere.reg.vtk

mris_convert ${subsurfdir}/rh.sphere.reg ${subsurfdir}/rh.sphere.reg.vtk



./BuildSphericalInterpolationMapping ${subsurfdir}/lh.sphere.reg.vtk ${refspheredir}/lh.sphere.ico5.reg.vtk ${subsurfdir}/lh.mapping.vtk

./BuildSphericalInterpolationMapping ${subsurfdir}/rh.sphere.reg.vtk ${refspheredir}/rh.sphere.ico5.reg.vtk ${subsurfdir}/rh.mapping.vtk

./ShapeInterpolation ${subsurfdir}/lh.white.vtk ${subsurfdir}/lh.mapping.vtk ${subsurfdir}/lh.white.corresponded.vtk

./ShapeInterpolation ${subsurfdir}/rh.white.vtk ${subsurfdir}/rh.mapping.vtk ${subsurfdir}/rh.white.corresponded.vtk

./ShapeInterpolation ${subsurfdir}/lh.midthickness.vtk ${subsurfdir}/lh.mapping.vtk ${subsurfdir}/lh.midthickness.corresponded.vtk

./ShapeInterpolation ${subsurfdir}/rh.midthickness.vtk ${subsurfdir}/rh.mapping.vtk ${subsurfdir}/rh.midthickness.corresponded.vtk
