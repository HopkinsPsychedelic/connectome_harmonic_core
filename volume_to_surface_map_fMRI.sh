#!/bin/bash


subsurfdir=$1

funcvolpath=$2

lhsavepath=$3

rhsavepath=$4










/home/neuro/repo/workbench-2/bin_rh_linux64/wb_command -volume-to-surface-mapping ${funcvolpath} ${subsurfdir}/lh.midthickness.corresponded.gii ${lhsavepath} -ribbon-constrained ${subsurfdir}/lh.white.corresponded.gii ${subsurfdir}/lh.pial.corresponded.gii



/home/neuro/repo/workbench-2/bin_rh_linux64/wb_command -volume-to-surface-mapping ${funcvolpath} ${subsurfdir}/rh.midthickness.corresponded.gii ${rhsavepath} -ribbon-constrained ${subsurfdir}/rh.white.corresponded.gii ${subsurfdir}/rh.pial.corresponded.gii


