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
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
#  THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# A stand-alone project that consumes an installed SELF the way a downstream
# application does. It is not part of SELF's own build; configure it against an
# install tree with
#
#   cmake -S test/external_consumer -B build-consumer -DCMAKE_PREFIX_PATH=<SELF install prefix>
#
# Nothing here adds SELF's include directory, definitions or dependencies by
# hand. If the program compiles, links and passes, SELF::self carried them.

cmake_minimum_required(VERSION 3.21)
project(SELFExternalConsumer LANGUAGES Fortran C)

find_package(SELF CONFIG REQUIRED)
message(STATUS "SELF ${SELF_VERSION} from ${SELF_DIR} : double precision ${SELF_DOUBLE_PRECISION}, GPU backend '${SELF_GPU_BACKEND}'")

add_executable(self_external_consumer self_external_consumer.f90)
set_source_files_properties(self_external_consumer.f90 PROPERTIES Fortran_PREPROCESS ON)
target_link_libraries(self_external_consumer PRIVATE SELF::self)

enable_testing()
add_test(NAME self_external_consumer COMMAND self_external_consumer)
