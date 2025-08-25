if [ -z ${CONDA_BUILD+x} ]; then
    source /opt/conda/conda-bld/roadies_1744817160316/work/build_env_setup.sh
fi
#!/bin/bash

# Debugging: Print current directory and list its contents
echo "Current directory: $(pwd)"
echo "Contents:"
ls -al

mkdir -p $PREFIX/ROADIES

# Build sampling code
if [[ ! -d "workflow/scripts/sampling/build" ]]; then
    cd workflow/scripts/sampling
    mkdir -p build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX="${PREFIX}"
    make -j"${CPU_COUNT}"
    cd ../../../..
fi

# Debugging: Print current directory and list its contents before copying
echo "Current directory before copying ROADIES: $(pwd)"
echo "Contents before copying ROADIES:"
ls -al

# Copy the entire ROADIES directory to the PREFIX directory
cp -rf * ${PREFIX}/ROADIES

# Debugging: Verify the contents of the PREFIX directory
echo "Contents of PREFIX:"
ls -al $PREFIX
