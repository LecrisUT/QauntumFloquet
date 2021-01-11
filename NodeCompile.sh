#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
# Make sure the following environment variables are created:
# PROJDIR=/path/to/project/root
# SLURM_JOB_DIR=/path/to/local/job/directory (Non-standard SLURM variable)

# Default variable preparation
if [ -z "${PROJDIR}" ]; then
  echo "PROJDIR environment variable was not set. Set it to the QuanFloq repository root accessible by the remote/node."
  exit 1
fi
if [ -z "${TARGET_DESTINATION}" ]; then
  TARGET_DESTINATION="${HOME}/${SLURM_JOB_PARTITION}"
fi
# Main Preparation
if [ ! -d ${SLURM_JOB_DIR} ]; then
  mkdir -p "${SLURM_JOB_DIR}"
fi
cd "${SLURM_JOB_DIR}" || exit 1
cp -R "${PROJDIR}/src" \
  "${PROJDIR}/test" \
  "${PROJDIR}/"*.cmake \
  "${PROJDIR}/CMakeLists.txt" \
  ./
mkdir build
# Main execution
cd build
cmake -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${TARGET_DESTINATION}" \
  -DBUILD_SHARED_LIBS=true \
  -DBUILD_VIRTUAL=False \
  -DMODULARITY=1 \
  -DDATA_TYPES=double,float,complex \
  -G "CodeBlocks - Unix Makefiles" ../
cmake --build ./ --target install -- -j "${SLURM_CPUS_ON_NODE}"
