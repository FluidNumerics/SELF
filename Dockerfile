FROM gcr.io/self-fluids/self-dep:latest AS devel

COPY . /tmp

RUN mkdir -p /tmp/extern

# FEQParse
RUN git clone https://github.com/FluidNumerics/feq-parse.git /tmp/extern/feq-parse && \
    mkdir -p /tmp/extern/feq-parse/build && \
    cd /tmp/extern/feq-parse/build && \
    cmake -DCMAKE_INSTALL_PREFIX="/apps/self" /tmp/extern/feq-parse && \
    make && make install

# JSON-Fortran
RUN git clone https://github.com/jacobwilliams/json-fortran.git /tmp/extern/json-fortran && \
    mkdir -p /tmp/extern/json-fortran/build && \
    cd /tmp/extern/json-fortran/build && \
    cmake -DSKIP_DOC_GEN=True -DCMAKE_INSTALL_PREFIX="/apps/self" /tmp/extern/json-fortran && \
    make && make install

# FLAP
RUN git clone --recurse-submodules https://github.com/szaghi/FLAP.git /tmp/extern/FLAP && \
    mkdir -p /tmp/extern/FLAP/build && \
    cd /tmp/extern/FLAP/build && \
    cmake -DCMAKE_INSTALL_PREFIX="/apps/self" /tmp/extern/FLAP && \
    make && make install

ENV HIP_PLATFORM=nvcc \
    HIP_COMPILER=/usr/local/cuda/bin/nvcc

RUN cd /tmp && \
    make install -f /tmp/self.make


FROM gcr.io/self-fluids/self-dep:latest

COPY --from=devel /apps/self /apps/self
ENV LD_LIBRARY_PATH=/usr/local/cuda/lib64:/apps/self/lib:$LD_LIBRARY_PATH \
    PATH=/usr/local/cuda/bin:/apps/self/bin:$PATH \
