#!/bin/bash
set -e

# For the results in the paper, every benchmark was run 50 times
ITERATIONS=50

CROSS_VARIANTS="rsdp rsdpg"
CROSS_CATEGORY="1 3 5"
CROSS_OPT="fast balanced small"

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
make PLATFORM=nucleo-l4r5zi -j $(nproc --all)

for variant in $CROSS_VARIANTS
do
    for category in $CROSS_CATEGORY
    do
        for opt in $CROSS_OPT
        do
            ./benchmarks.py --platform nucleo-l4r5zi --uart /dev/ttyACM1 -i $ITERATIONS --nohashing cross-$variant-$category-$opt
        done
    done
done
./convert_benchmarks.py md
