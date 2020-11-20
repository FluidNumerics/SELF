FROM gcr.io/self-fluids/rocm:latest AS devel

# FEQParse
RUN git clone https://github.com/FluidNumerics/feq-parse.git /tmp/feq-parse &&\
    mkdir -p /tmp/feq-parse/build && \
    cd /tmp/feq-parse/build && \
    cmake -DCMAKE_INSTALL_PREFIX=/apps/self /tmp/feq-parse && \
    make && make install && \
    rm -rf /tmp/feq-parse

# JSON-Fortran
RUN git clone https://github.com/jacobwilliams/json-fortran.git /tmp/json-fortran &&\
    mkdir -p /tmp/json-fortran/build && \
    cd /tmp/json-fortran/build && \
    cmake -DSKIP_DOC_GEN=True -DCMAKE_INSTALL_PREFIX=/apps/self /tmp/json-fortran && \
    make && make install &&\
    rm -rf /tmp/json-fortran

ENV LD_LIBRARY_PATH=/apps/self/lib:$LD_LIBRARY_PATH \
    PATH=/apps/self/bin:$PATH

RUN source /opt/rh/devtoolset-9/enable && \
    mkdir -p /self-build/build && \
    cd /self-build/build && \
    FC="/opt/rocm/bin/hipfc" \
    CC="/opt/rocm/bin/hipfc" \
    CXX="/opt/rocm/bin/hipfc" \
    FFLAGS="-v -DGPU -ffree-line-length-none" \
    CXXFLAGS="-v" \
    /usr/local/bin/cmake /self-build -DCMAKE_INSTALL_PREFIX=/usr/local/self &&\
    make && make install


FROM gcr.io/self-fluids/rocm:latest

COPY --from=devel /apps/self
