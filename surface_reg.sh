#!/bin/bash

subsurfdir = $1
refdir = $2




mris_expand -thickness ${subsurfdir}/lh.white 0.5 ${subsurfdir}/surf/lh.mid

mris_expand -thickness ${subsurfdir}/surf/rh.white 0.5 ${subsurfdir}/rh.mid

mris_convert ${subsurfdir}/lh.mid ${subsurfdir}/lh.mid.vtk

mris_convert ${subsurfdir}/rh.mid ${subsurfdir}/rh.mid.vtk

mris_convert ${subsurfdir}/lh.sphere ${subsurfdir}/lh.sphere.gii
mris_convert ${subsurfdir}/rh.sphere ${subsurfdir}/rh.sphere.gii

mris_convert ${refdir}/lh.sphere ${refdir}/lh.sphere.gii
mris_convert ${refdir}/rh.sphere ${refdir}/rh.sphere.gii

mris_convert ${refdir}/lh.sphere ${refdir}/lh.sphere.vtk
mris_convert ${refdir}/rh.sphere ${refdir}/rh.sphere.vtk

mri_convert ${subsurfdir}/lh.curv ${subsurfdir}/lh.curv.gii
mri_convert ${subsurfdir}/rh.curv ${subsurfdir}/rh.curv.gii

mri_convert ${subsurfdir}/lh.sulc ${subsurfdir}/lh.sulc.gii
mri_convert ${subsurfdir}/rh.sulc ${subsurfdir}/rh.sulc.gii

mri_convert ${refdir}/lh.sulc ${refdir}/lh.sulc.gii
mri_convert ${refdir}/rh.sulc ${refdir}/rh.sulc.gii

msm --inmesh=${subsurfdir}/lh.sphere.gii --refmesh=${refdir}/lh.sphere.gii --indata=${subsurfdir}/lh.sulc.gii --refdata=${refdir}/lh.sulc.gii -o ${subsurfdir}/lh.MSM.sphere.reg.surf.gii

msm --inmesh=${subsurfdir}/rh.sphere.gii --refmesh=${refdir}/rh.sphere.gii --indata=${subsurfdir}/rh.sulc.gii --refdata=${refdir}/rh.sulc.gii -o ${subsurfdir}/rh.MSM.sphere.reg.surf.gii

mris_convert ${subsurfdir}/lh.MSM.sphere.reg.surf.gii ${subsurfdir}/lh.MSM.sphere.reg.surf.vtk
mris_convert ${subsurfdir}/rh.MSM.sphere.reg.surf.gii ${subsurfdir}/rh.MSM.sphere.reg.surf.vtk

./BuildSphericalInterpolationMapping ${subsurfdir}/lh.MSM.sphere.reg.surf.vtk ${refdir}/lh.sphere.vtk ${subsurfdir}/lh.mapping.vtk

./BuildSphericalInterpolationMapping ${subsurfdir}/rh.MSM.sphere.reg.surf.vtk ${refdir}/rh.sphere.vtk ${subsurfdir}/lh.mapping.vtk

./ShapeInterpolation ${subsurfdir}/lh.white.vtk ${subsurfdir}/lh.mapping.vtk ${subsurfdir}/lh.white.fsaverage.corresponded.vtk

./ShapeInterpolation ${subsurfdir}/rh.white.vtk ${subsurfdir}/rh.mapping.vtk ${subsurfdir}/rh.white.fsaverage.corresponded.vtk
