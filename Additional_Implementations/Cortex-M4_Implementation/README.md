# Cortex-M4 implementation of CROSS as presented in the paper "Optimizing the Post Quantum Signature Scheme CROSS for Resource Constrained Devices" [2]

## Summary
This folder contains a highly configurable runtime memory optimized implementation of the Codes and Restricted Objects Signature Scheme (CROSS). Our build scripts by default generate four different implementations based on these optimizations and export them to the pqm4 framework.

Namely, these implementations are:
* `ref`: The reference implementation adjusted for pqm4. It contains an additional header file to export build-time defined options to pqm4's build process.
* `m4stack`: A runtime memory optimized implementation with all stack saving options of this work enabled.
* `m4speed`: A runtime optimized version for Cortex-M4. As some of the runtime optimizations increase the stack size the are not by default enabled.
* `m4opt`: A version with all stack and speed optimizations enabled.

One can tune these implementations by enabling or disabling any of the optimizations independently as discussed below. The implementations above are configured as discussed in the paper.

## Dependencies
On Debian-based distributions, the required dependencies can be installed as follows:
```
sudo apt install git cmake unifdef python3 python3-pip python3-venv
```
Though untested, on other Linux distributions or MacOS, please adjust the above command according to your platform, on Windows, running the benchmarks should be possible using [WSL](https://learn.microsoft.com/de-de/windows/wsl/install).

We originally used python3.12, `cmake` `3.28.3` and `unifdef` `2.12` for this work.
One furthermore needs a compiler for ARM Cortex-m4, throughout the paper we used `arm-gcc (arm-none-eabi)` version [13.2](https://developer.arm.com/downloads/-/arm-gnu-toolchain-downloads).
This artifact was furthermore tested with version 14.2 and 15.2 of `arm-gcc`. Please note however that not all other implementations in `pqm4` compile with these newer versions.
One furthermore needs either `stlink` or `openocd` (depending on the target platform) to program the binaries, the [README of pqm4](https://github.com/mupq/pqm4) contains the necessary instructions.

## Required Hardware
The benchmarking results for this work were taken on the Nucleo-l4r5zi board as outlined in the documentation of pqm4 [1]. Additional results were furthermore obtained on the STM32F4Discovery board as also outlined in the documentation of pqm4.

## Quick start
In case one just wants to generate the implementations used throughout this paper and run them in the pqm4 framework on their default device we provided a script. This assumes that all required dependencies are installed and that the necessary benchmarking device is connected as outlined in [1]. Please note that the script is currently configured to run each benchmark 50 times as done for the paper which takes substantial time (several hours). This iteration count can be configured at the top of the script to e.g. 1 iteration for testing purposes though this slightly increases variability on the results. Please note that the serial interface for communication with the controller is likely to be found under a different/name path on your platform and adjust the variable at the top of the script accordingly. Reproducing the results obtained on the STM32F4 is also possible by adjusting the platform and serial path at the top of the script.

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

### Switching different optimizations on/off
To switch a specific optimization on/off one only needs to set it 0 or 1 in `Cortex-M4_Implementation/CMakeLists.txt`.

### Exporting the implementations to pqm4
When running `cmake ..` from `Cortex-M4_Implementation/build` CMake exports the configured implementations to `Cortex-M4_Implementation/build/pqm4`. Please follow the instructions from pqm4 from here onward or have a look at `run.sh`.

## Repositories
This implementation, including possible future updates, can also be found at the [Official CROSS Github Repository](https://github.com/CROSS-signature/CROSS-implementation).
The configured implementations for pqm4 as used for this work are also available at [pqm4](https://github.com/mupq/pqm4), possibly as open merge request from [here](https://github.com/joschupp/pqm4).

## References
[1]: [pqm4](https://github.com/mupq/pqm4)
[2]: [eprint version of the paper](https://eprint.iacr.org/2025/1928)
