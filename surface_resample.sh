#!/bin/bash

subsurfdir=$1

refspheredir=$2






mris_convert ${subsurfdir}/lh.midthickness ${subsurfdir}/lh.midthickness.vtk

mris_convert ${subsurfdir}/rh.midthickness ${subsurfdir}/rh.midthickness.vtk

mris_convert ${subsurfdir}/lh.white ${subsurfdir}/lh.white.vtk

mris_convert ${subsurfdir}/rh.white ${subsurfdir}/rh.white.vtk

mris_convert ${subsurfdir}/lh.pial ${subsurfdir}/lh.pial.vtk

mris_convert ${subsurfdir}/rh.pial ${subsurfdir}/rh.pial.vtk

mris_convert ${subsurfdir}/lh.sphere.reg ${subsurfdir}/lh.sphere.reg.vtk

mris_convert ${subsurfdir}/rh.sphere.reg ${subsurfdir}/rh.sphere.reg.vtk



/home/neuro/repo/BuildSphericalInterpolationMapping ${subsurfdir}/lh.sphere.reg.vtk ${refspheredir}/lh.sphere.ico5.reg.vtk ${subsurfdir}/lh.mapping.vtk

/home/neuro/repo/BuildSphericalInterpolationMapping ${subsurfdir}/rh.sphere.reg.vtk ${refspheredir}/rh.sphere.ico5.reg.vtk ${subsurfdir}/rh.mapping.vtk

/home/neuro/repo/ShapeInterpolation ${subsurfdir}/lh.white.vtk ${subsurfdir}/lh.mapping.vtk ${subsurfdir}/lh.white.corresponded.vtk

/home/neuro/repo/ShapeInterpolation ${subsurfdir}/rh.white.vtk ${subsurfdir}/rh.mapping.vtk ${subsurfdir}/rh.white.corresponded.vtk

/home/neuro/repo/ShapeInterpolation ${subsurfdir}/lh.midthickness.vtk ${subsurfdir}/lh.mapping.vtk ${subsurfdir}/lh.midthickness.corresponded.vtk

/home/neuro/repo/ShapeInterpolation ${subsurfdir}/rh.midthickness.vtk ${subsurfdir}/rh.mapping.vtk ${subsurfdir}/rh.midthickness.corresponded.vtk

/home/neuro/repo/ShapeInterpolation ${subsurfdir}/lh.pial.vtk ${subsurfdir}/lh.mapping.vtk ${subsurfdir}/lh.pial.corresponded.vtk

/home/neuro/repo/ShapeInterpolation ${subsurfdir}/rh.pial.vtk ${subsurfdir}/rh.mapping.vtk ${subsurfdir}/rh.pial.corresponded.vtk


mris_convert ${subsurfdir}/lh.midthickness.corresponded.vtk ${subsurfdir}/lh.midthickness.corresponded.gii

mris_convert ${subsurfdir}/rh.midthickness.corresponded.vtk ${subsurfdir}/rh.midthickness.corresponded.gii

mris_convert ${subsurfdir}/lh.white.corresponded.vtk ${subsurfdir}/lh.white.corresponded.gii

mris_convert ${subsurfdir}/rh.white.corresponded.vtk ${subsurfdir}/rh.white.corresponded.gii

mris_convert ${subsurfdir}/lh.pial.corresponded.vtk ${subsurfdir}/lh.pial.corresponded.gii

mris_convert ${subsurfdir}/rh.pial.corresponded.vtk ${subsurfdir}/rh.pial.corresponded.gii