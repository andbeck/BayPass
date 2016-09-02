#!/bin/bash

boomShackalaka () {
	
	efile_suffix=${1}
	
	echo ${efile_suffix}
	
	g_baypass \
	-npop 8 \
	-gfile ALLELEFILE \
	-efile ENVFILE_${eFile_suffix} \
	-poolsizefile SAMPLEFILE \
	—audxmodel \
	-d0yij 20 \
	-outprefix ana${eFile_suffix}_
}
export -f boomShackalaka

parallel -j1 boomShackalaka ::: BASE BLat BlatTemp BlatTempPH




