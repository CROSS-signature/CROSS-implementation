# Cortex-M4 implementation of CROSS as presented in the paper "Optimizing the Post Quantum Signature Scheme CROSS for Resource Constrained Devices" [2]

## Summary
This folder contains a highly configurable runtime memory optimized implementation of the Codes and Restricted Objects Signature Scheme (CROSS). Our build scripts by default generate four different implementations based on these optimizations and export them to the pqm4 framework.

Namely, these implementations are:
* `ref`: The reference implementation adjusted for pqm4. It contains an additional header file to export build-time defined otpions to pqm4's build process.
* `m4stack`: A runtime memory optimized implementation with all stack svaing options of this work enabled.
* `m4speed`: A runtime optimized version for Cortex-M4. As some of the runtime optimizations increase the stack size the are not by default enabled.
* `m4opt`: A version with all stack and speed optimizations enabled.

One can tune these implementations by enabling or disabling any of the optimizations independently as discussed below. The implementations above are configured as discussed in the paper.

## Dependencies
On Debian-based distributions, the requried dependencies can be installed as follows (please adjust for other systems accordingly):
```
sudo apt install git cmake unifdef python3 python3-pip python3-venv
```
We originally used python3.12 for this work.
One furthermore needs a compiler for ARM Cortex-m4, throughout this work we used `arm-gcc` version 13.2.

## Required Hardware
The benchmarking results for this work were taken on the Nucleo-l4r5zi board as outlined in the documentation of pqm4 [1].

## Quick start
In case one just wants to generate the implementations used throughout this paper and run them in the pqm4 framework on their default device we provided a script. This assumes that all required dependencies are installed and that the necessary benchmarking device is connected as outlined in [1]. Please note that the script is currently configured to run each benchmark 50 times as done for the paper which takes substantial time (several hours). This iteration count can be configured at the top of the script to e.g. 1 iteration for testing purposes though this slightly increases variability on the results.

```bash
./run.sh
```

# Structure of this folder
* `m4ref`: Contains files adapted from the original cross reference for compatibility with the pqm4 framework (inclusion of randombytes as dependency)
* `m4speed`: Contains the source files specific for the runtime optimizations on Cortex-M4
* `m4stack`: Contains the source files specific for the stack optimizations on Cortex-M4
* `pqm4`: Contains the files necessary for interaction with pqm4, namely a header to export build definitions as well as a patch for the skiplist. This is *not* the pqm4-folder into which the pqm4 repo is cloned and into which the implementations get exported.
* `build`: To be created by the user for cmake as outlined above.
* `CMakeLists.txt`: CMake configuration for exporting the implementations to pqm4
* `test_all_mem_opts.py`: A script to test the functional correctness of all available optimization configurations that don't rely on ARM assembly on x86.
* `run.sh`: A simple script exporting the implementation configurations used throughout this paper to pqm4 and benchmarking them on their default target device.

## Detailed Usage

### General
This implementation can be configured using different preprocessor directives. As the CROSS reference implementation uses cmake as build-framework, we employ the same framework to enable/disable these directives.
Based on the configured options in `Cortex-M4_Implementation/CMakeLists.txt` we then export these implementations with the configured options to `pqm4`.
As `pqm4` ships its own build framework we opted to use `unifdef` to reduce the exported source files to the parts enabled by the chosen optimizations.
As the parameter set, optimization corner and variant of the CROSS reference implementation is usually also configured using CMake, we export an additional `build_defs.h` to make these settings available in pqm4.

### Switching different optimiations on/off
To switch a specific optimization on/off one only needs to set it 0 or 1 in `Cortex-M4_Implementation/CMakeLists.txt`.

### Exporting the implementations to pqm4
When running `cmake ..` from `Cortex-M4_Implementation/build` CMake exports the configured implementations to `Cortex-M4_Implementation/build/pqm4`. Please follow the instructions from pqm4 from here onwards or have a look at `run.sh`.

## Repositories
This implementation, including possible future updates, can also be found at the [Offical CROSS Github Repository](https://github.com/CROSS-signature/CROSS-implementation).
The configured implementations for pqm4 as used for this work are also available at [pqm4](https://github.com/mupq/pqm4), possibly as open merge request from [here](https://github.com/joschupp/pqm4).

## References
\[1\]: [pqm4](https://github.com/mupq/pqm4)
\[2\]: [eprint version of the paper](https://eprint.iacr.org/2025/1928)
