#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
# 
#  Maintainers : support@fluidnumerics.com
#  Official Repository : https://github.com/FluidNumerics/self/
# 
#  Copyright © 2024 Fluid Numerics LLC
# 
#  Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 
#  1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 
#  2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in
#     the documentation and/or other materials provided with the distribution.
# 
#  3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from
#     this software without specific prior written permission.
# 
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
#  THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

CMAKE_MINIMUM_REQUIRED(VERSION 3.21)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
set(make_env SELF_HOME=${CMAKE_BINARY_DIR} SELF_COMPILER=${SELF_COMPILER})

function (add_fortran_tests)

    foreach (TFILE ${ARGN})
        string (REGEX REPLACE ".f90$" "" TEST_NAME ${TFILE})
        add_executable (${TEST_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/${TFILE})
        add_dependencies(${TEST_NAME} self)
        target_link_libraries(${TEST_NAME} self)
        target_include_directories(${TEST_NAME} PUBLIC ${CMAKE_BINARY_DIR}/include)
        add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
	set_tests_properties(${TEST_NAME} PROPERTIES LABELS "serial")
        install(TARGETS ${TEST_NAME} DESTINATION test)
    endforeach ()

endfunction ()

function (add_mpi_fortran_tests)

    foreach (TFILE ${ARGN})
        string (REGEX REPLACE ".f90$" "" TEST_NAME ${TFILE})
        add_executable (${TEST_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/${TFILE})
        add_dependencies(${TEST_NAME} self)
        target_link_libraries(${TEST_NAME} self)
        target_include_directories(${TEST_NAME} PUBLIC ${CMAKE_BINARY_DIR}/include)
        add_test(NAME ${TEST_NAME} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS} ${SELF_MPIEXEC_OPTIONS} ${TEST_NAME})
	set_tests_properties(${TEST_NAME} PROPERTIES LABELS "parallel")

        install(TARGETS ${TEST_NAME} DESTINATION test)
    endforeach ()

endfunction ()

add_fortran_tests (
    "scalargridinterp_1d_constant.f90"
    "scalarderivative_1d_constant.f90"
    "scalarboundaryinterp_1d_constant.f90"
    "scalargridinterp_2d_constant.f90"
    "scalargradient_2d_constant.f90"
    "scalarboundaryinterp_2d_constant.f90"
    "vectordivergence_2d_constant.f90"
    "vectordivergence_2d_linear.f90"
    "vectorboundaryinterp_2d_constant.f90"
    "scalargridinterp_3d_constant.f90"
    "scalargradient_3d_constant.f90"
    "scalarboundaryinterp_3d_constant.f90"
    "vectordivergence_3d_constant.f90"
    "vectorgridinterp_2d_constant.f90"
    "vectorgridinterp_3d_constant.f90"
    "vectorboundaryinterp_3d_constant.f90"
    "mesh2d_setup.f90"
    "mesh2d_uniformstructured.f90"
    "mesh2d_uniform_refine.f90"
    "mesh3d_uniform_refine.f90"
    "solution_transfer_2d.f90"
    "solution_transfer_3d.f90"
    "solution_transfer_2d_device.f90"
    "solution_transfer_3d_device.f90"
    "solution_transfer_2d_device_offset.f90"
    "solution_transfer_3d_device_offset.f90"
    "quadtree_balance_2d.f90"
    "octree_balance_3d.f90"
    "adaptive_mortar_2d.f90"
    "adaptive_mortar_3d.f90"
    "transfer_plan_2d.f90"
    "transfer_plan_3d.f90"
    "transfer_plan_2d_guard_badleaf.f90"
    "transfer_plan_3d_guard_badleaf.f90"
    "transfer_plan_2d_guard_snapshot.f90"
    "transfer_plan_3d_guard_snapshot.f90"
    "transfer_plan_2d_guard_unbuilt.f90"
    "transfer_plan_3d_guard_unbuilt.f90"
    "transfer_plan_2d_guard_range.f90"
    "transfer_plan_3d_guard_range.f90"
    "transfer_plan_2d_window.f90"
    "transfer_plan_3d_window.f90"
    "transfer_plan_2d_guard_window.f90"
    "transfer_plan_3d_guard_window.f90"
    "transfer_plan_2d_guard_windowbounds.f90"
    "transfer_plan_3d_guard_windowbounds.f90"
    "transfer_plan_2d_guard_windowfamily.f90"
    "transfer_plan_3d_guard_windowfamily.f90"
    "amr_controller_2d_guard_decomp.f90"
    "amr_controller_3d_guard_decomp.f90"
    "lineareuler2d_amr_soundwave.f90"
    "lineareuler3d_amr_soundwave.f90"
    "lineareuler2d_amr_coarsen_wake.f90"
    "dgmodel2d_regrid_guard_uninit.f90"
    "dgmodel3d_regrid_guard_uninit.f90"
    "data_2d_device_memory.f90"
    "geometry_2d_reuse.f90"
    "geometry_3d_reuse.f90"
    "scalar_2d_guard_resize_nvar.f90"
    "scalar_3d_guard_resize_nvar.f90"
    "vector_2d_guard_resize_nvar.f90"
    "vector_3d_guard_resize_nvar.f90"
    "tensor_2d_guard_resize_nvar.f90"
    "tensor_3d_guard_resize_nvar.f90"
    "dgmodel2d_guard_transfer_unstaged.f90"
    "dgmodel3d_guard_transfer_unstaged.f90"
    "dgmodel2d_guard_window_and_stage.f90"
    "dgmodel3d_guard_window_and_stage.f90"
    "dgmodel2d_guard_download_nowindow.f90"
    "dgmodel3d_guard_download_nowindow.f90"
    "dgmodel2d_guard_download_badbuffer.f90"
    "dgmodel3d_guard_download_badbuffer.f90"
    "dgmodel3d_guard_transfer_reapply.f90"
    "amr_controller_2d_guard_maxlevel.f90"
    "amr_controller_3d_guard_maxlevel.f90"
    "amr_controller_2d_guard_wrongmesh.f90"
    "amr_controller_3d_guard_wrongmesh.f90"
    "refinement_indicator_2d_guard_thresholds.f90"
    "refinement_indicator_3d_guard_thresholds.f90"
    "refinement_indicator_2d_guard_setthresholds.f90"
    "refinement_indicator_3d_guard_setthresholds.f90"
    "refinement_indicator_2d_guard_degree.f90"
    "refinement_indicator_3d_guard_degree.f90"
    "refinement_indicator_2d_guard_nelem.f90"
    "refinement_indicator_3d_guard_nelem.f90"
    "refinement_indicator_2d_guard_ivar.f90"
    "refinement_indicator_3d_guard_ivar.f90"
    "refinement_indicator_2d_guard_relfloor.f90"
    "refinement_indicator_3d_guard_relfloor.f90"
    "refinement_indicator_2d_guard_energyscale.f90"
    "refinement_indicator_3d_guard_energyscale.f90"
    "refinement_indicator_2d_guard_energyweights.f90"
    "refinement_indicator_3d_guard_energyweights.f90"
    "refinement_indicator_2d_guard_negweight.f90"
    "refinement_indicator_3d_guard_negweight.f90"
    "refinement_indicator_2d_guard_emptyweights.f90"
    "refinement_indicator_3d_guard_emptyweights.f90"
    "refinement_indicator_2d_guard_weightcount.f90"
    "refinement_indicator_3d_guard_weightcount.f90"
    "refinement_indicator_2d_guard_gatesize.f90"
    "refinement_indicator_3d_guard_gatesize.f90"
    "refinement_indicator_2d_guard_significantfloor.f90"
    "refinement_indicator_3d_guard_significantfloor.f90"
    "adaptive_mesh_2d_guard_unbalanced.f90"
    "adaptive_mesh_3d_guard_unbalanced.f90"
    "meshrefine_2d_guard_multirank.f90"
    "meshrefine_3d_guard_multirank.f90"
    "quadtree_2d_guard_multirank.f90"
    "octree_3d_guard_multirank.f90"
    "mesh2d_read_ism.f90"
    "mesh2d_read_ismmm.f90"
    "mesh3d_read_ism.f90"
    "mesh3d_read_ismmm.f90"
    "mesh3d_read_manymaterials.f90"
    "mesh3d_read_rotatedstack.f90"
    "mesh3d_read_guard_missingfile.f90"
    "mesh3d_read_guard_badconnectivity.f90"
    "mesh3d_setup.f90"
    "mesh3d_uniformstructured.f90"
    "mappedscalarderivative_1d_constant.f90"
    "mappedscalarbrderivative_1d_constant.f90"
    "mappedscalardgderivative_1d_constant.f90"
    "mappedscalargradient_2d_constant.f90"
    "mappedscalargradient_2d_linear.f90"
    "mappedscalarbrgradient_2d_constant.f90"
    "mappedscalarbrgradient_2d_constant_structuredmesh.f90"
    "mappedscalarbrgradient_2d_linear.f90"
    "mappedscalarbrgradient_2d_linear_structuredmesh.f90"
    "mappedvectordivergence_2d_constant.f90"
    "mappedvectordivergence_2d_gausslobatto_constant.f90"
    "mappedvectordivergence_2d_linear_structuredmesh.f90"
    "mappedvectordgdivergence_2d_constant.f90"
    "mappedvectordgdivergence_2d_linear.f90"
    "mappedvectordgdivergence_2d_linear_structuredmesh.f90"
    "mortarprojection_identity.f90"
    "refinement_indicator_2d_spectraldecay.f90"
    "refinement_indicator_3d_spectraldecay.f90"
    "refinement_indicator_2d_relative_floor.f90"
    "refinement_indicator_3d_relative_floor.f90"
    "refinement_indicator_2d_energy_hysteresis.f90"
    "refinement_indicator_3d_energy_hysteresis.f90"
    "mappedscalarmortarexchange_2d_linear.f90"
    "mappedscalarmortarexchange_3d_linear.f90"
    "mappedvectordgdivergence_2d_mortar.f90"
    "mappedvectordgdivergence_3d_mortar.f90"
    "mappedscalargradient_3d_constant.f90"
    "mappedscalargradient_3d_linear.f90"
    "mappedscalarbrgradient_3d_constant.f90"
    "mappedscalarbrgradient_3d_constant_structuredmesh.f90"
    "mappedscalarbrgradient_3d_linear.f90"
    "mappedscalarbrgradient_3d_linear_n8.f90"
    "mappedscalarbrgradient_3d_linear_structuredmesh.f90"
    "mappedvectordivergence_3d_constant.f90"
    "mappedvectordivergence_3d_linear.f90"
    "mappedvectordivergence_3d_linear_n8.f90"
    "mappedvectordivergence_3d_gausslobatto_constant.f90"
    "mappedvectordgdivergence_3d_constant.f90"
    "mappedvectordgdivergence_3d_linear.f90"
    "mappedvectordgdivergence_3d_linear_n8.f90"
    "mappedvectordgdivergence_3d_linear_n12.f90"
    "mappedvectordgdivergence_3d_linear_structuredmesh.f90"
    "mappedvectordgdivergence_3d_linear_sideexchange.f90"
    "nulldgmodel1d.f90"
    "nulldgmodel2d.f90"
    "nulldgmodel2d_prescribed.f90"
    "nulldgmodel3d.f90"
    "nulldgmodel3d_prescribed.f90"
    "advection_diffusion_1d_euler.f90"
    "advection_diffusion_1d_euler_pickup.f90"
    "advection_diffusion_1d_rk2.f90"
    "advection_diffusion_1d_rk3.f90"
    "advection_diffusion_1d_rk4.f90"
    "burgers1d_constant.f90"
    "burgers1d_nonormalflow.f90"
    "burgers1d_prescribed.f90"
    "burgers1d_radiation.f90"
    "advection_diffusion_2d_euler.f90"
    "advection_diffusion_2d_rk2.f90"
    "advection_diffusion_2d_rk3.f90"
    "advection_diffusion_2d_rk3_pickup.f90"
    "advection_diffusion_2d_rk4.f90"
    "advection_diffusion_2d_nostress.f90"
    "advection_diffusion_2d_mortar.f90"
    "advection_diffusion_3d_euler.f90"
    "advection_diffusion_3d_rk2.f90"
    "advection_diffusion_3d_rk3.f90"
    "advection_diffusion_3d_rk4.f90"
    "advection_diffusion_3d_nostress.f90"
    "linear_shallow_water_2d_constant.f90"
    "linear_shallow_water_2d_nonormalflow.f90"
    "linear_shallow_water_2d_radiation.f90"
    "lineareuler2d_pml_absorption.f90"
    "lineareuler2d_mortar_soundwave.f90"
    "lineareuler3d_soundwave.f90"
    "lineareuler3d_mortar_soundwave.f90"
    "twopointvectordivergence_2d_constant.f90"
    "twopointvectordivergence_2d_linear.f90"
    "twopointvectordivergence_2d_rotational.f90"
    "twopointvectordivergence_3d_constant.f90"
    "twopointvectordivergence_3d_linear.f90"
    "mappedtwopointvectordivergence_2d_constant.f90"
    "mappedtwopointvectordivergence_2d_linear.f90"
    "mappedtwopointvectordivergence_3d_constant.f90"
    "mappedtwopointvectordivergence_3d_linear.f90"
    "ec_advection_2d_rk3.f90"
    "ec_advection_2d_entropy_conservation.f90"
    "ec_advection_2d_mortar_guard.f90"
    "ec_advection_3d_mortar_guard.f90"
    "ec_advection_3d_rk3.f90"
    "ec_advection_3d_entropy_conservation.f90"
    "esatmo3d_motionless.f90"
    "esatmo3d_hydrostaticbalance.f90"
    "points_locate_2d.f90"
    "points_locate_2d_fine.f90"
    "points_locate_2d_outside.f90"
    "points_locate_3d.f90"
    "points_eval_cache_2d.f90"
    "points_eval_cache_3d.f90"
    "points_eval_gpu_2d.f90"
    "points_eval_gpu_3d.f90"
    "points_dirac_delta_2d.f90"
    "points_dirac_delta_3d.f90"
    "points_dirac_delta_2d_ndim_mismatch.f90"
    "points_dirac_delta_3d_ndim_mismatch.f90"
    "points_dirac_delta_2d_nvar_mismatch.f90"
    "points_dirac_delta_3d_nvar_mismatch.f90"
    "domaindecomposition_two_meshes.f90"
    )

# Mark the four DiracDelta guard tests as expected-to-fail: each program
# intentionally triggers a stop 1 inside the guard, and a non-zero exit
# code is the success criterion.
set_tests_properties(
    points_dirac_delta_2d_ndim_mismatch
    points_dirac_delta_3d_ndim_mismatch
    points_dirac_delta_2d_nvar_mismatch
    points_dirac_delta_3d_nvar_mismatch
    ec_advection_2d_mortar_guard
    ec_advection_3d_mortar_guard
    refinement_indicator_2d_guard_thresholds
    refinement_indicator_3d_guard_thresholds
    refinement_indicator_2d_guard_setthresholds
    refinement_indicator_3d_guard_setthresholds
    refinement_indicator_2d_guard_degree
    refinement_indicator_3d_guard_degree
    refinement_indicator_2d_guard_nelem
    refinement_indicator_3d_guard_nelem
    refinement_indicator_2d_guard_ivar
    refinement_indicator_3d_guard_ivar
    refinement_indicator_2d_guard_relfloor
    refinement_indicator_3d_guard_relfloor
    refinement_indicator_2d_guard_energyscale
    refinement_indicator_3d_guard_energyscale
    refinement_indicator_2d_guard_energyweights
    refinement_indicator_3d_guard_energyweights
    refinement_indicator_2d_guard_negweight
    refinement_indicator_3d_guard_negweight
    refinement_indicator_2d_guard_emptyweights
    refinement_indicator_3d_guard_emptyweights
    refinement_indicator_2d_guard_weightcount
    refinement_indicator_3d_guard_weightcount
    refinement_indicator_2d_guard_gatesize
    refinement_indicator_3d_guard_gatesize
    refinement_indicator_2d_guard_significantfloor
    refinement_indicator_3d_guard_significantfloor
    adaptive_mesh_2d_guard_unbalanced
    adaptive_mesh_3d_guard_unbalanced
    meshrefine_2d_guard_multirank
    meshrefine_3d_guard_multirank
    quadtree_2d_guard_multirank
    octree_3d_guard_multirank
    transfer_plan_2d_guard_badleaf
    transfer_plan_3d_guard_badleaf
    transfer_plan_2d_guard_snapshot
    transfer_plan_3d_guard_snapshot
    transfer_plan_2d_guard_unbuilt
    transfer_plan_3d_guard_unbuilt
    transfer_plan_2d_guard_range
    transfer_plan_3d_guard_range
    transfer_plan_2d_guard_window
    transfer_plan_3d_guard_window
    transfer_plan_2d_guard_windowbounds
    transfer_plan_3d_guard_windowbounds
    transfer_plan_2d_guard_windowfamily
    transfer_plan_3d_guard_windowfamily
    amr_controller_2d_guard_decomp
    amr_controller_3d_guard_decomp
    dgmodel2d_regrid_guard_uninit
    dgmodel3d_regrid_guard_uninit
    amr_controller_2d_guard_maxlevel
    amr_controller_3d_guard_maxlevel
    amr_controller_2d_guard_wrongmesh
    amr_controller_3d_guard_wrongmesh
    scalar_2d_guard_resize_nvar
    scalar_3d_guard_resize_nvar
    vector_2d_guard_resize_nvar
    vector_3d_guard_resize_nvar
    tensor_2d_guard_resize_nvar
    tensor_3d_guard_resize_nvar
    dgmodel2d_guard_transfer_unstaged
    dgmodel3d_guard_transfer_unstaged
    dgmodel2d_guard_window_and_stage
    dgmodel3d_guard_window_and_stage
    dgmodel2d_guard_download_nowindow
    dgmodel3d_guard_download_nowindow
    dgmodel2d_guard_download_badbuffer
    dgmodel3d_guard_download_badbuffer
    dgmodel3d_guard_transfer_reapply
    mesh3d_read_guard_missingfile
    mesh3d_read_guard_badconnectivity
    PROPERTIES WILL_FAIL TRUE)

add_mpi_fortran_tests( "mappedvectordgdivergence_2d_linear_mpi.f90"
                       "mappedvectordgdivergence_2d_linear_sideexchange_mpi.f90"
                       "mappedvectordgdivergence_2d_linear_structuredmesh_mpi.f90"
                       "mappedscalarbrgradient_2d_linear_mpi.f90"
                       "mappedvectordgdivergence_3d_linear_mpi.f90"
                       "mappedvectordgdivergence_3d_linear_sideexchange_mpi.f90"
                       "mappedvectordgdivergence_3d_linear_structuredmesh_mpi.f90"
                       "mappedscalarbrgradient_3d_linear_mpi.f90"
                       "mappedscalarmortarexchange_2d_linear_mpi.f90"
                       "mappedscalarmortarexchange_3d_linear_mpi.f90"
                       "mappedvectordgdivergence_2d_mortar_mpi.f90"
                       "mappedvectordgdivergence_3d_mortar_mpi.f90"
                       "lineareuler2d_mortar_soundwave_mpi.f90"
                       "lineareuler3d_mortar_soundwave_mpi.f90"
                       "lineareuler2d_amr_soundwave_mpi.f90"
                       "lineareuler3d_amr_soundwave_mpi.f90"
                       "amr_migrate_2d_window_mpi.f90"
    "amr_migrate_2d_device_apply_mpi.f90"
                       "amr_migrate_3d_window_mpi.f90"
    "amr_migrate_3d_device_apply_mpi.f90"
                       "geometry_2d_reuse_mpi.f90"
                       "advection_diffusion_2d_mortar_mpi.f90"
                       "advection_diffusion_2d_rk3_mpi.f90"
                       "advection_diffusion_2d_rk3_pickup_mpi.f90"
                       "advection_diffusion_3d_rk3_mpi.f90"
                       "advection_diffusion_3d_rk3_pickup_mpi.f90"
                       "structuredmesh_external_comm_2d.f90"
                       "structuredmesh_external_comm_3d.f90"
                       "structuredmesh_subcomm_2d.f90"
                       "read_hopr_external_comm_2d.f90" )

# The Stage 6c geometry diagnostics (SELF_AMR_GEOM_*) are env-gated, so ordinary runs never enter
# them. Re-run existing AMR tests with each switch set: this exercises the verification and bypass
# paths for real rather than leaving them as untested code, and it costs no new sources because the
# executables are already built above.
#
#  - GEOM_VERIFY   cross-checks the incrementally assembled geometry against full regeneration
#  - GEOM_FULL     bypasses the incremental assembly (must reproduce the same physics)
#  - GEOM_NO_REUSE regenerates every element (the reuse predicate never fires)
#
# lineareuler2d_amr_soundwave asserts per-epoch conservation and non-increasing entropy, so it is a
# real check under each switch and not merely a coverage exercise.
add_test(NAME geometry_2d_reuse_mpi_verify
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS}
                 ${SELF_MPIEXEC_OPTIONS} geometry_2d_reuse_mpi)
set_tests_properties(geometry_2d_reuse_mpi_verify PROPERTIES
                     LABELS "parallel" ENVIRONMENT "SELF_AMR_GEOM_VERIFY=1")

add_test(NAME lineareuler2d_amr_soundwave_geom_full COMMAND lineareuler2d_amr_soundwave)
set_tests_properties(lineareuler2d_amr_soundwave_geom_full PROPERTIES
                     LABELS "serial" ENVIRONMENT "SELF_AMR_GEOM_FULL=1")

add_test(NAME lineareuler2d_amr_soundwave_geom_noreuse COMMAND lineareuler2d_amr_soundwave)
set_tests_properties(lineareuler2d_amr_soundwave_geom_noreuse PROPERTIES
                     LABELS "serial" ENVIRONMENT "SELF_AMR_GEOM_NO_REUSE=1")

# 3-D equivalents: the same switches gate the 3-D controller's incremental geometry path, and
# lineareuler3d_amr_soundwave carries the same conservation/entropy assertions.
add_test(NAME lineareuler3d_amr_soundwave_geom_verify COMMAND lineareuler3d_amr_soundwave)
set_tests_properties(lineareuler3d_amr_soundwave_geom_verify PROPERTIES
                     LABELS "serial" ENVIRONMENT "SELF_AMR_GEOM_VERIFY=1")

add_test(NAME lineareuler3d_amr_soundwave_geom_full COMMAND lineareuler3d_amr_soundwave)
set_tests_properties(lineareuler3d_amr_soundwave_geom_full PROPERTIES
                     LABELS "serial" ENVIRONMENT "SELF_AMR_GEOM_FULL=1")

add_test(NAME lineareuler3d_amr_soundwave_geom_noreuse COMMAND lineareuler3d_amr_soundwave)
set_tests_properties(lineareuler3d_amr_soundwave_geom_noreuse PROPERTIES
                     LABELS "serial" ENVIRONMENT "SELF_AMR_GEOM_NO_REUSE=1")

# Multi-rank solution migration (issue #167). The default path is now the point-to-point window
# exchange, so the MPI AMR tests above already cover it. These re-runs cover the two things the
# default run cannot show, again with no new sources:
#
#  - MIGRATE_VERIFY  runs BOTH migrations in one process and compares the received window against
#                    the allgathered global old field bit for bit. This is the assertion that the
#                    point-to-point routing is right, not merely self-consistent.
#  - MIGRATE_GATHER  keeps the retained v1 allgather path exercised, so the documented fallback
#                    cannot rot.
add_test(NAME lineareuler2d_amr_soundwave_mpi_migrate_verify
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS}
                 ${SELF_MPIEXEC_OPTIONS} lineareuler2d_amr_soundwave_mpi)
set_tests_properties(lineareuler2d_amr_soundwave_mpi_migrate_verify PROPERTIES
                     LABELS "parallel" ENVIRONMENT "SELF_AMR_MIGRATE_VERIFY=1")

add_test(NAME lineareuler3d_amr_soundwave_mpi_migrate_verify
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS}
                 ${SELF_MPIEXEC_OPTIONS} lineareuler3d_amr_soundwave_mpi)
set_tests_properties(lineareuler3d_amr_soundwave_mpi_migrate_verify PROPERTIES
                     LABELS "parallel" ENVIRONMENT "SELF_AMR_MIGRATE_VERIFY=1")

add_test(NAME lineareuler2d_amr_soundwave_mpi_migrate_gather
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS}
                 ${SELF_MPIEXEC_OPTIONS} lineareuler2d_amr_soundwave_mpi)
set_tests_properties(lineareuler2d_amr_soundwave_mpi_migrate_gather PROPERTIES
                     LABELS "parallel" ENVIRONMENT "SELF_AMR_MIGRATE_GATHER=1")

add_test(NAME lineareuler3d_amr_soundwave_mpi_migrate_gather
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS}
                 ${SELF_MPIEXEC_OPTIONS} lineareuler3d_amr_soundwave_mpi)
set_tests_properties(lineareuler3d_amr_soundwave_mpi_migrate_gather PROPERTIES
                     LABELS "parallel" ENVIRONMENT "SELF_AMR_MIGRATE_GATHER=1")

# Re-runs of the same sources for the windowed transfer (#172). TRANSFER_VERIFY applies the plan
# BOTH ways from the same migrated window and compares - to a TOLERANCE, and only because the device
# compiler contracts the kernel's multiply-accumulates into FMAs. It is deliberately a DIFFERENT
# switch from MIGRATE_VERIFY, which asserts the window ARRIVED intact and is exact; merging them
# would drag the strict one down to the loose bar. Registered in both dimensions because
# VerifyWindowedApply is a per-dimension routine, and a twin covering the same logic is not a
# substitute for running it.
add_test(NAME lineareuler2d_amr_soundwave_mpi_transfer_verify
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS}
                 ${SELF_MPIEXEC_OPTIONS} lineareuler2d_amr_soundwave_mpi)
set_tests_properties(lineareuler2d_amr_soundwave_mpi_transfer_verify PROPERTIES
                     LABELS "parallel" ENVIRONMENT "SELF_AMR_TRANSFER_VERIFY=1")

add_test(NAME lineareuler3d_amr_soundwave_mpi_transfer_verify
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS}
                 ${SELF_MPIEXEC_OPTIONS} lineareuler3d_amr_soundwave_mpi)
set_tests_properties(lineareuler3d_amr_soundwave_mpi_transfer_verify PROPERTIES
                     LABELS "parallel" ENVIRONMENT "SELF_AMR_TRANSFER_VERIFY=1")

# TRANSFER_HOST forces the portable host migration and windowed apply even on a GPU build. On a CPU
# build it selects the same code the default path already takes, so what it covers here is the branch
# selection itself - but that branch is what routes a GPU build back to the portable path, and
# leaving it unexecuted anywhere means a typo in it would only surface on hardware. 2-D only; the
# 3-D twin is the same three lines.
add_test(NAME lineareuler2d_amr_soundwave_mpi_transfer_host
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS}
                 ${SELF_MPIEXEC_OPTIONS} lineareuler2d_amr_soundwave_mpi)
set_tests_properties(lineareuler2d_amr_soundwave_mpi_transfer_host PROPERTIES
                     LABELS "parallel" ENVIRONMENT "SELF_AMR_TRANSFER_HOST=1")

# Both switches together, which is not redundant with either alone: it is the only thing that
# reaches the transferHost branch of the WINDOW DOWNLOAD. On a GPU build that branch reads the
# window out of host memory while the default reads it off the device, so a bug there would
# silently corrupt what MIGRATE_VERIFY and TRANSFER_VERIFY then compare against - the diagnostics
# would be lying rather than failing.
add_test(NAME lineareuler2d_amr_soundwave_mpi_transfer_host_verify
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS}
                 ${SELF_MPIEXEC_OPTIONS} lineareuler2d_amr_soundwave_mpi)
set_tests_properties(lineareuler2d_amr_soundwave_mpi_transfer_host_verify PROPERTIES
                     LABELS "parallel"
                     ENVIRONMENT "SELF_AMR_TRANSFER_HOST=1;SELF_AMR_TRANSFER_VERIFY=1")

# Two ranks cannot distinguish gather-then-slice from a point-to-point routing bug: with one peer
# the send and receive runs are forced, and no rank's window can span more than one peer. These
# entries re-launch the SAME executables on SELF_MPIEXEC_NUMPROCS_MANY ranks, where a window can
# straddle several peers and a rank can end up with nothing to own.
add_test(NAME lineareuler2d_amr_soundwave_mpi_np_many
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS_MANY}
                 ${SELF_MPIEXEC_OPTIONS} ${SELF_MPIEXEC_OPTIONS_MANY}
                 lineareuler2d_amr_soundwave_mpi)
set_tests_properties(lineareuler2d_amr_soundwave_mpi_np_many PROPERTIES
                     LABELS "parallel" PROCESSORS ${SELF_MPIEXEC_NUMPROCS_MANY})

add_test(NAME lineareuler2d_amr_soundwave_mpi_np_many_migrate_verify
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS_MANY}
                 ${SELF_MPIEXEC_OPTIONS} ${SELF_MPIEXEC_OPTIONS_MANY}
                 lineareuler2d_amr_soundwave_mpi)
set_tests_properties(lineareuler2d_amr_soundwave_mpi_np_many_migrate_verify PROPERTIES
                     LABELS "parallel" PROCESSORS ${SELF_MPIEXEC_NUMPROCS_MANY}
                     ENVIRONMENT "SELF_AMR_MIGRATE_VERIFY=1")

add_test(NAME amr_migrate_2d_window_mpi_np_many
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS_MANY}
                 ${SELF_MPIEXEC_OPTIONS} ${SELF_MPIEXEC_OPTIONS_MANY}
                 amr_migrate_2d_window_mpi)
set_tests_properties(amr_migrate_2d_window_mpi_np_many PROPERTIES
                     LABELS "parallel" PROCESSORS ${SELF_MPIEXEC_NUMPROCS_MANY})

add_test(NAME amr_migrate_3d_window_mpi_np_many
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS_MANY}
                 ${SELF_MPIEXEC_OPTIONS} ${SELF_MPIEXEC_OPTIONS_MANY}
                 amr_migrate_3d_window_mpi)
set_tests_properties(amr_migrate_3d_window_mpi_np_many PROPERTIES
                     LABELS "parallel" PROCESSORS ${SELF_MPIEXEC_NUMPROCS_MANY})

# The model-level device migration and windowed apply, at more than two ranks. Two ranks force every
# send/receive run and cannot produce a window that straddles more than one peer, so the 4-rank
# entries are where the routing is actually distinguishable from gather-then-slice.
add_test(NAME amr_migrate_2d_device_apply_mpi_np_many
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS_MANY}
                 ${SELF_MPIEXEC_OPTIONS} ${SELF_MPIEXEC_OPTIONS_MANY}
                 amr_migrate_2d_device_apply_mpi)
set_tests_properties(amr_migrate_2d_device_apply_mpi_np_many PROPERTIES
                     LABELS "parallel" PROCESSORS ${SELF_MPIEXEC_NUMPROCS_MANY})

add_test(NAME amr_migrate_3d_device_apply_mpi_np_many
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS_MANY}
                 ${SELF_MPIEXEC_OPTIONS} ${SELF_MPIEXEC_OPTIONS_MANY}
                 amr_migrate_3d_device_apply_mpi)
set_tests_properties(amr_migrate_3d_device_apply_mpi_np_many PROPERTIES
                     LABELS "parallel" PROCESSORS ${SELF_MPIEXEC_NUMPROCS_MANY})

add_test(NAME lineareuler3d_amr_soundwave_mpi_np_many_migrate_verify
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS_MANY}
                 ${SELF_MPIEXEC_OPTIONS} ${SELF_MPIEXEC_OPTIONS_MANY}
                 lineareuler3d_amr_soundwave_mpi)
set_tests_properties(lineareuler3d_amr_soundwave_mpi_np_many_migrate_verify PROPERTIES
                     LABELS "parallel" PROCESSORS ${SELF_MPIEXEC_NUMPROCS_MANY}
                     ENVIRONMENT "SELF_AMR_MIGRATE_VERIFY=1")
