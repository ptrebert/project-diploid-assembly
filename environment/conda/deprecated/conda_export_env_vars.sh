#!/usr/bin/env bash

REPOSITORY_PREFIX="/home/pebert/work/code/github/project-diploid-assembly"

CONDA_PREFIX="/TL/epigenetics2/work/pebert/conda/envs/dipassm"

cd ${CONDA_PREFIX}
mkdir -p ./etc/conda/activate.d
mkdir -p ./etc/conda/deactivate.d
#touch ./etc/conda/activate.d/env_vars.sh
#touch ./etc/conda/deactivate.d/env_vars.sh

cp -f ${REPOSITORY_PREFIX}/environment/conda/activate/env_vars.sh ./etc/conda/activate.d/env_vars.sh
cp -f ${REPOSITORY_PREFIX}/environment/conda/deactivate/env_vars.sh ./etc/conda/deactivate.d/env_vars.sh