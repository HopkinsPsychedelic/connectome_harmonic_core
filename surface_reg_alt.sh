#!/bin/bash

subsurfdir=$1
refdir=$2
lrefsphere=$3
rrefsphere=$3

lrefsulc=$4
rrefsulc=$5







mris_expand -thickness ${subsurfdir}/lh.white 0.5 ${subsurfdir}/lh.mid

mris_expand -thickness ${subsurfdir}/surf/rh.white 0.5 ${subsurfdir}/rh.mid

mris_convert ${subsurfdir}/lh.mid ${subsurfdir}/lh.mid.vtk

mris_convert ${subsurfdir}/rh.mid ${subsurfdir}/rh.mid.vtk

mris_convert ${refdir}/${lrefsphere} ${refdir}/${lrefsphere}.gii
mris_convert ${refdir}/${rrefsphere} ${refdir}/${rrefsphere}.gii

mris_convert ${subsurfdir}/lh.sphere ${subsurfdir}/lh.sphere.gii
mris_convert ${subsurfdir}/rh.sphere ${subsurfdir}/rh.sphere.gii

mris_convert ${refdir}/${lrefsphere} ${refdir}/${lrefsphere}.vtk
mris_convert ${refdir}/${rrefsphere} ${refdir}/${rrefsphere}.vtk

mri_convert ${subsurfdir}/lh.curv ${subsurfdir}/lh.curv.gii
mri_convert ${subsurfdir}/rh.curv ${subsurfdir}/rh.curv.gii

mri_convert ${subsurfdir}/lh.sulc ${subsurfdir}/lh.sulc.gii
mri_convert ${subsurfdir}/rh.sulc ${subsurfdir}/rh.sulc.gii

mri_convert ${refdir}/${lrefsulc} ${refdir}/${lrefsulc}.gii
mri_convert ${refdir}/${rrefsulc} ${refdir}/${rrefsulc}.gii


msm --inmesh=${subsurfdir}/lh.sphere.gii --refmesh=${refdir}/${lrefsphere}.gii --indata=${subsurfdir}/lh.sulc.gii --refdata=${refdir}/${lrefsulc}.gii -o ${subsurfdir}/lh.MSM.sphere.reg.surf.gii

msm --inmesh=${subsurfdir}/rh.sphere.gii --refmesh=${refdir}/${rrefsphere}.gii --indata=${subsurfdir}/rh.sulc.gii --refdata=${refdir}/${rrefsulc}.gii -o ${subsurfdir}/rh.MSM.sphere.reg.surf.gii

mris_convert ${subsurfdir}/lh.MSM.sphere.reg.surf.gii ${subsurfdir}/lh.MSM.sphere.reg.surf.vtk
mris_convert ${subsurfdir}/rh.MSM.sphere.reg.surf.gii ${subsurfdir}/rh.MSM.sphere.reg.surf.vtk

./BuildSphericalInterpolationMapping ${subsurfdir}/lh.MSM.sphere.reg.surf.vtk ${refdir}/${lrefsphere}.vtk ${subsurfdir}/lh.mapping.vtk

./BuildSphericalInterpolationMapping ${subsurfdir}/rh.MSM.sphere.reg.surf.vtk ${refdir}/${rrefsphere}.vtk ${subsurfdir}/lh.mapping.vtk

./ShapeInterpolation ${subsurfdir}/lh.white.vtk ${subsurfdir}/lh.mapping.vtk ${subsurfdir}/lh.white.fsaverage.corresponded.vtk

./ShapeInterpolation ${subsurfdir}/rh.white.vtk ${subsurfdir}/rh.mapping.vtk ${subsurfdir}/rh.white.fsaverage.corresponded.vtk
