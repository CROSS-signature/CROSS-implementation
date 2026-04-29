cmake_minimum_required(VERSION 3.7)

project(CROSS C)
set(CMAKE_C_STANDARD 11)

# Reference code dir, same across all CROSS implementations
set(REFERENCE_CODE_DIR ../../Reference_Implementation)
# Code dir w. adaptions for pqm4
set(PQM4_CODE_DIR ./pqm4)
set(M4_STACK_OPT_CODE_DIR ./m4stack)
set(M4_SPEED_OPT_CODE_DIR ./m4speed)
set(M4_REF_CODE_DIR ./m4ref)
set(M4_TESTING_CODE_DIR ./testing)

set(PQM4_PATH pqm4/crypto_sign/)
set(MUPQ_PATH pqm4/mupq/crypto_sign/)
# set(PQM4_PATH ../../../../../../pqm4/crypto_sign/)
# set(MUPQ_PATH ../../../../../../pqm4/mupq/crypto_sign/)

set(TEST_MODE 0)

if (CUSTOM_OPT_TARGET)
    set(TEST_MODE 1)
endif()

if (!OPT_STRING)
    set(OPT_STRING "")
endif()

# Sources and Headers from reference, optimized further down
set(BASE_DIR ${REFERENCE_CODE_DIR})
set(SPEC_HEADERS_ref
        ${M4_STACK_OPT_CODE_DIR}/include/api.h
        ${BASE_DIR}/include/sha3.h
        ${BASE_DIR}/include/restr_arith.h
        ${BASE_DIR}/include/fp_arith.h
        ${BASE_DIR}/include/seedtree.h
        ${BASE_DIR}/include/merkle_tree.h
        ${BASE_DIR}/include/pack_unpack.h
        ${M4_REF_CODE_DIR}/include/csprng_hash.h
)
set(SPEC_SOURCES_ref
        ${BASE_DIR}/lib/merkle.c
        ${BASE_DIR}/lib/seedtree.c
        ${M4_STACK_OPT_CODE_DIR}/lib/sign.c
        ${BASE_DIR}/lib/csprng_hash.c
        ${BASE_DIR}/lib/pack_unpack.c
        ${M4_REF_CODE_DIR}/lib/CROSS.c
)

set(SPEC_SOURCES_ref_testing
        ${BASE_DIR}/lib/keccakf1600.c
        ${BASE_DIR}/lib/fips202.c
)

# Sources and Headers from m4stack opt.
set(BASE_DIR ${M4_STACK_OPT_CODE_DIR})
set(SPEC_HEADERS_m4stack
        ${BASE_DIR}/include/api.h
        ${BASE_DIR}/include/sha3.h
        ${BASE_DIR}/include/csprng_hash.h
        ${BASE_DIR}/include/restr_arith.h
        ${BASE_DIR}/include/fp_arith.h
        ${BASE_DIR}/include/seedtree.h
        ${BASE_DIR}/include/merkle_tree.h
        ${BASE_DIR}/include/m4_defs.h
        ${BASE_DIR}/include/pack_unpack.h
)
set(SPEC_SOURCES_m4stack
        ${BASE_DIR}/lib/merkle.c
        ${BASE_DIR}/lib/seedtree.c
        ${BASE_DIR}/lib/CROSS.c
        ${BASE_DIR}/lib/sign.c
        ${BASE_DIR}/lib/csprng_hash.c
        ${BASE_DIR}/lib/pack_unpack.c
)

set(BASE_DIR ${M4_SPEED_OPT_CODE_DIR})
set(SPEC_HEADERS_m4speed
        ${SPEC_HEADERS_m4stack}
        ${BASE_DIR}/include/fp_arith.S
)
set(SPEC_SOURCES_m4speed ${SPEC_SOURCES_m4stack})

set(SPEC_HEADERS_m4opt ${SPEC_HEADERS_m4speed})
set(SPEC_SOURCES_m4opt ${SPEC_SOURCES_m4speed})

set(COMMON_DIR ${REFERENCE_CODE_DIR})
set(HEADERS
    ${COMMON_DIR}/include/CROSS.h
    ${COMMON_DIR}/include/parameters.h
)

set(TESTING_HEADERS
    ${M4_TESTING_CODE_DIR}/include/randombytes.h)

# Memory optimization flags sign operation
set(MEM_OPT_SIGN_INC_CMT_1 1)
set(MEM_OPT_SIGN_RECOMP_E_V_U 1)
set(MEM_OPT_SIGN_RECOMP_Y 1)
set(MEM_OPT_SIGN_INC_MTREE 1)
set(MEM_OPT_SIGN_INC_STREE 1)
set(MEM_OPT_SIGN_BS_SEEDS_FROM_TREE 1)

# Memory optimization flags verify operation
set(MEM_OPT_VERIFY_INC_CMT_1 1)
set(MEM_OPT_VERIFY_INC_Y 1)
set(MEM_OPT_VERIFY_INC_MTREE 1)
set(MEM_OPT_VERIFY_INC_STREE 1)
set(MEM_OPT_VERIFY_BS_SEEDS_FROM_TREE 1)

# Memory optimization for all operations
set(MEM_OPT_HASH_TO_MTREE 1)
set(MEM_OPT_IN_PLACE_SAMP 1)
# set(MEM_OPT_INC_SAMPLING 1)
# set(MEM_OPT_STREE_IN_CCM 1)
# set(MEM_OPT_MTREE_IN_CCM 1)

# Time optimization flags
# set(TIME_OPT_SIMD 1)
# Table Look-Up for exponentation, keep this disabled unless you know exactly what you are doing!
set(TIME_OPT_TBL_LOOKUP 0)

set(OPTS MEM_OPT_SIGN_INC_CMT_1 MEM_OPT_SIGN_RECOMP_E_V_U MEM_OPT_SIGN_RECOMP_Y MEM_OPT_SIGN_INC_MTREE MEM_OPT_SIGN_INC_STREE MEM_OPT_SIGN_BS_SEEDS_FROM_TREE MEM_OPT_VERIFY_INC_CMT_1 MEM_OPT_VERIFY_INC_Y MEM_OPT_VERIFY_INC_MTREE MEM_OPT_VERIFY_INC_STREE MEM_OPT_VERIFY_BS_SEEDS_FROM_TREE MEM_OPT_HASH_TO_MTREE MEM_OPT_IN_PLACE_SAMP MEM_OPT_INC_SAMPLING MEM_OPT_STREE_IN_CCM MEM_OPT_MTREE_IN_CCM TIME_OPT_TBL_LOOKUP)

if(TEST_MODE EQUAL 0 AND NOT GEN_KAT EQUAL 1)
    message("Cloning the official pqm4 repo to build/pqm4")
    execute_process(COMMAND bash -c "[ \"$(ls -A pqm4)\" ] && echo \"The pqm4 folder already exists, assuming correct checkout\" || git clone --recursive https://github.com/mupq/pqm4")
    message("Remove existing (Round 1) CROSS implementations from pqm4")
    execute_process(COMMAND bash -c "rm -r pqm4/crypto_sign/cross*")
    execute_process(COMMAND bash -c "rm -r pqm4/mupq/crypto_sign/cross*")
    set(IMPLEMENTATIONS ref m4stack m4speed m4opt)
    # set(IMPLEMENTATIONS m4speed)
    foreach(impl ${IMPLEMENTATIONS})
        set(UNIFDEF " ")
		if(impl STREQUAL m4speed OR impl STREQUAL m4opt)
            set(TIME_OPT_SIMD 1)
            set(UNIFDEF "${UNIFDEF} -DTIME_OPT_SIMD")
        else()
            unset(TIME_OPT_SIMD)
            set(UNIFDEF "${UNIFDEF} -UTIME_OPT_SIMD")
        endif()
        if(impl STREQUAL ref)
            set(COPY_PATH ${MUPQ_PATH})
        else()
            set(COPY_PATH ${PQM4_PATH})
        endif()
		if(impl STREQUAL m4speed)
			foreach(OPT ${OPTS})
				set(UNIFDEF "${UNIFDEF} -U${OPT}")
			endforeach(OPT)
		else()
			foreach(OPT ${OPTS})
				if(${OPT} EQUAL 1)
					set(UNIFDEF "${UNIFDEF} -D${OPT}")
				else()
					set(UNIFDEF "${UNIFDEF} -U${OPT}")
				endif()
			endforeach(OPT)
		endif()
        message(STATUS "${UNIFDEF}")
        set(COPIED 0)
        foreach(category RANGE 1 5 2)
            set(RSDP_VARIANTS RSDP RSDPG)
            foreach(RSDP_VARIANT ${RSDP_VARIANTS})
				string(TOLOWER ${RSDP_VARIANT} rsdp_lower)
                set(PARAM_TARGETS SIG_SIZE BALANCED SPEED)
                foreach(optimiz_target ${PARAM_TARGETS})
                    unset(NO_TREES)
                    if(optimiz_target STREQUAL SPEED)
                        set(NO_TREES 1)
                    endif()
					if(optimiz_target STREQUAL SPEED)
						set(opt_name fast)
					elseif(optimiz_target STREQUAL BALANCED)
						set(opt_name balanced)
					elseif(optimiz_target STREQUAL SIG_SIZE)
						set(opt_name small)
					endif()
                    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory cross-${rsdp_lower}-${category}-${opt_name}/${impl} WORKING_DIRECTORY ${COPY_PATH})
                    configure_file(${PQM4_CODE_DIR}/include/build_defs_unifdef.h ${COPY_PATH}cross-${rsdp_lower}-${category}-${opt_name}/${impl}/build_defs.h)
                    SET(SOURCES_CPY ${SPEC_SOURCES_${impl}})
                    foreach(source_file ${SOURCES_CPY})
                        get_filename_component(filename ${source_file} NAME)
                        if(COPIED EQUAL 0)
                            configure_file(${source_file} ${COPY_PATH}cross-${rsdp_lower}-${category}-${opt_name}/${impl}/${filename} COPYONLY)
                            execute_process(COMMAND bash -c "unifdef ${UNIFDEF} -o ${COPY_PATH}cross-${rsdp_lower}-${category}-${opt_name}/${impl}/${filename} ${COPY_PATH}cross-${rsdp_lower}-${category}-${opt_name}/${impl}/${filename}")
                        else()
                            execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${TARGET_PATH}/${filename} ${COPY_PATH}cross-${rsdp_lower}-${category}-${opt_name}/${impl}/${filename})
                        endif()
                    endforeach(source_file)
                    SET(HEADERS_CPY ${HEADERS} ${SPEC_HEADERS_${impl}})
                    foreach(header_file ${HEADERS_CPY})
                        get_filename_component(filename ${header_file} NAME)
                        if(COPIED EQUAL 0)
                            configure_file(${header_file} ${COPY_PATH}cross-${rsdp_lower}-${category}-${opt_name}/${impl}/${filename} COPYONLY)
                            execute_process(COMMAND bash -c "unifdef ${UNIFDEF} -o ${COPY_PATH}cross-${rsdp_lower}-${category}-${opt_name}/${impl}/${filename} ${COPY_PATH}cross-${rsdp_lower}-${category}-${opt_name}/${impl}/${filename}")
                        else()
                            execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${TARGET_PATH}/${filename} ${COPY_PATH}cross-${rsdp_lower}-${category}-${opt_name}/${impl}/${filename})
                        endif()
                    endforeach(header_file)
                    if(COPIED EQUAL 0)
                        set(TARGET_PATH ../../cross-${rsdp_lower}-${category}-${opt_name}/${impl})
                    endif()
                    set(COPIED 1)
                endforeach(optimiz_target)
            endforeach(RSDP_VARIANT)
        endforeach(category)
    endforeach(impl)
elseif(NOT GEN_KAT EQUAL 1)
    message("Testing Mode, testing all optimization combinations...")
    message(${CUSTOM_OPT_TARGET})
    foreach(category RANGE 1 5 2)
        set(RSDP_VARIANTS RSDP RSDPG)
        foreach(RSDP_VARIANT ${RSDP_VARIANTS})
            set(PARAM_TARGETS SIG_SIZE BALANCED SPEED)
            foreach(optimiz_target ${PARAM_TARGETS})
                set(impl m4stack)
                unset(NO_TREES)
                if(optimiz_target STREQUAL SPEED)
                    set(NO_TREES 1)
                endif()
                if(optimiz_target STREQUAL SPEED)
                    set(OMIT_SEED_TREE "-DNO_TREES=1")
                else()
                    set(OMIT_SEED_TREE "")
                endif()
                SET(SOURCES_CPY ${SPEC_SOURCES_${impl}} ${SPEC_SOURCES_ref_testing})
                SET(HEADERS_CPY ${HEADERS} ${SPEC_HEADERS_${impl}} ${TESTING_HEADERS})
                set(TARGET_BINARY_NAME CROSS_test_cat_${category}_${RSDP_VARIANT}_${optimiz_target}_${OPT_STRING})
                add_executable(${TARGET_BINARY_NAME} ${HEADERS_CPY} ${SPEC_SOURCES_${impl}} ${SPEC_SOURCES_ref_testing} ../Benchmarking/include/arith_unit_tests.h
                               ../Benchmarking/lib/CROSS_test.c)
                target_include_directories(${TARGET_BINARY_NAME} PRIVATE
                                           ${BASE_DIR}/include
                                           ${M4_STACK_OPT_CODE_DIR}/include
                                           ${COMMON_DIR}/include
                                           ${M4_TESTING_CODE_DIR}/include
                                           ../Benchmarking/include)
                target_link_libraries(${TARGET_BINARY_NAME} m ${SANITIZE} ${KECCAK_EXTERNAL_LIB})
                set_target_properties(${TARGET_BINARY_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./bin)
                set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY
                             COMPILE_FLAGS "${OMIT_SEED_TREE} ${CUSTOM_OPT_TARGET} -DCATEGORY_${category} -D${optimiz_target}=1 -D${RSDP_VARIANT}=1 ${KECCAK_EXTERNAL_ENABLE} ")
            endforeach(optimiz_target)
        endforeach(RSDP_VARIANT)
    endforeach(category)
else()
    message("Testing Mode, building for KAT generation...")
    # message(${CUSTOM_OPT_TARGET})
    foreach(category RANGE 1 5 2)
        set(RSDP_VARIANTS RSDP RSDPG)
        foreach(RSDP_VARIANT ${RSDP_VARIANTS})
            set(PARAM_TARGETS SIG_SIZE BALANCED SPEED)
            foreach(optimiz_target ${PARAM_TARGETS})
                set(impl m4stack)
                unset(NO_TREES)
                if(optimiz_target STREQUAL SPEED)
                    set(NO_TREES 1)
                endif()
                if(optimiz_target STREQUAL SPEED)
                    set(OMIT_SEED_TREE "-DNO_TREES=1")
                else()
                    set(OMIT_SEED_TREE "")
                endif()
                SET(SOURCES_CPY ${SPEC_SOURCES_${impl}} ${SPEC_SOURCES_ref_testing})
                SET(HEADERS_CPY ${HEADERS} ${SPEC_HEADERS_${impl}} ${TESTING_HEADERS})
                set(TARGET_BINARY_NAME CROSS_KATgen_cat_${category}_${RSDP_VARIANT}_${optimiz_target}_${OPT_STRING})
                add_executable(${TARGET_BINARY_NAME} ${HEADERS_CPY} ${SPEC_SOURCES_${impl}} ${SPEC_SOURCES_ref_testing}
                               ../KAT_Generation/lib/rng.c ../KAT_Generation/lib/PQCgenKAT_sign.c)
                target_include_directories(${TARGET_BINARY_NAME} PRIVATE
                                           ${BASE_DIR}/include
                                           ${M4_STACK_OPT_CODE_DIR}/include
                                           ${COMMON_DIR}/include
                                           ${M4_TESTING_CODE_DIR}/include
                                           ../KAT_Generation/include)
                target_link_libraries(${TARGET_BINARY_NAME} m ssl crypto ${SANITIZE} ${KECCAK_EXTERNAL_LIB})
                set_target_properties(${TARGET_BINARY_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./bin)
                set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY
                             COMPILE_FLAGS "${OMIT_SEED_TREE} ${CUSTOM_OPT_TARGET} -DCATEGORY_${category} -D${optimiz_target}=1 -D${RSDP_VARIANT}=1 ${KECCAK_EXTERNAL_ENABLE} ")
            endforeach(optimiz_target)
        endforeach(RSDP_VARIANT)
    endforeach(category)
endif()
