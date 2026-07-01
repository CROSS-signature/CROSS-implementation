#!/bin/bash
set -e

# For the results in the paper, every benchmark was run 50 times
ITERATIONS=50

# The results presented in the paper were obtained on either the nucleo-l4r5zi or on the stm32f4discovery
PLATFORM=nucleo-l4r5zi
# PLATFORM=stm32f4discovery

# The serial interface used by pqm4 for the benchmarks likely needs to be adapted to your platform,
# on Ubuntu/Debian it is usually /dev/ttyACM* for the nucleo-l4r5zi 
# and some /dev/ttyUSB* for the external adapter required for the stm32f4discovery
SERIAL_PATH=/dev/ttyACM1

CROSS_VARIANTS="rsdp rsdpg"
CROSS_CATEGORY="1 3 5"
CROSS_OPT="fast balanced small"
CROSS_IMPLS="m4stack m4speed m4opt"

cd Cortex-M4_Implementation
mkdir -p build
cd build
cmake ..
cd pqm4
python3 -m venv pqm4_venv
. pqm4_venv/bin/activate
pip install -r requirements.txt
if git apply --check ../../pqm4/skiplist.patch; then
    git apply ../../pqm4/skiplist.patch
else
    echo "The pqm4 skiplist patch can not be applied, either this has already been done or needs to be done manually if pqm4 has been updated in between, please check."
fi
if git apply --check ../../pqm4/keccakf1600.patch; then
    git apply ../../pqm4/keccakf1600.patch
else
    echo "The pqm4 keccakf1600 patch can not be applied, either this has already been done or needs to be done manually if pqm4 has been updated in between, please check."
fi

# Build only the implementations relevent for this paper as some of the other schemes have issues with later compiler versions
for variant in $CROSS_VARIANTS
do
    for category in $CROSS_CATEGORY
    do
        for opt in $CROSS_OPT
        do
            make -j $(nproc --all) PLATFORM=$PLATFORM IMPLEMENTATION_PATH=mupq/crypto_sign/cross-$variant-$category-$opt/ref
            for impl in $CROSS_IMPLS
            do
                make -j $(nproc --all) PLATFORM=$PLATFORM IMPLEMENTATION_PATH=crypto_sign/cross-$variant-$category-$opt/$impl
            done
        done
    done
done

for variant in $CROSS_VARIANTS
do
    for category in $CROSS_CATEGORY
    do
        for opt in $CROSS_OPT
        do
            ./benchmarks.py --platform $PLATFORM --uart $SERIAL_PATH -i $ITERATIONS --nohashing cross-$variant-$category-$opt
        done
    done
done
./convert_benchmarks.py md
