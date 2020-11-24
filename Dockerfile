FROM gcr.io/self-fluids/self-dep:latest AS devel

COPY . /tmp

# FEQParse
RUN mkdir -p /tmp/extern/feq-parse/build && \
    cd /tmp/extern/feq-parse/build && \
    cmake -DCMAKE_INSTALL_PREFIX="/apps/feq-parse" /tmp/extern/feq-parse && \
    make && make install

# JSON-Fortran
RUN mkdir -p /tmp/extern/json-fortran/build && \
    cd /tmp/extern/json-fortran/build && \
    cmake -DSKIP_DOC_GEN=True -DCMAKE_INSTALL_PREFIX="/apps/json-fortran" /tmp/extern/json-fortran && \
    make && make install

# FLAP
RUN mkdir -p /tmp/extern/FLAP/build && \
    cd /tmp/extern/FLAP/build && \
    cmake -DCMAKE_INSTALL_PREFIX="/apps/flap" /tmp/extern/FLAP && \
    make && make install

ENV LD_LIBRARY_PATH=/apps/feq-parse/lib:/apps/flap/lib:$LD_LIBRARY_PATH

RUN mkdir -p /tmp/build && \
    cd /tmp/build && \
    FC="/opt/hipfort/bin/hipfc" \
    CC="/opt/hipfort/bin/hipfc" \
    CXX="/opt/hipfort/bin/hipfc" \
    FFLAGS="-v -DGPU -ffree-line-length-none" \
    CXXFLAGS="-v" \
    cmake -DCMAKE_INSTALL_PREFIX="/apps/self" /tmp &&\
    export HIP_PLATFORM=nvcc && export HIP_COMPILER=/usr/local/cuda/bin/nvcc && make VERBOSE=1 && make install


FROM gcr.io/self-fluids/self-dep:latest

COPY --from=devel /apps /apps
