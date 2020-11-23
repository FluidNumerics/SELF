FROM gcr.io/self-fluids/self-dep:latest AS devel

COPY . /tmp

# FEQParse
RUN mkdir -p /tmp/dependencies/feq-parse/build && \
    cd /tmp/dependencies/feq-parse/build && \
    cmake -DCMAKE_INSTALL_PREFIX="/apps/self" /tmp/dependencies/feq-parse && \
    make && make install

# JSON-Fortran
RUN mkdir -p /tmp/dependencies/json-fortran/build && \
    cd /tmp/dependencies/json-fortran/build && \
    cmake -DSKIP_DOC_GEN=True -DCMAKE_INSTALL_PREFIX="/apps/self" /tmp/dependencies/json-fortran && \
    make && make install

ENV LD_LIBRARY_PATH=/apps/self/lib:$LD_LIBRARY_PATH \
    PATH=/apps/self/bin:$PATH

RUN mkdir -p /tmp/build && \
    cd /tmp/build && \
    FC="/usr/local/bin/hipfc" \
    CC="/usr/local/bin/hipfc" \
    CXX="/usr/local/bin/hipfc" \
    FFLAGS="-v -DGPU -ffree-line-length-none" \
    CXXFLAGS="-v" \
    cmake -DCMAKE_INSTALL_PREFIX="/apps/self" /tmp &&\
    make && make install


FROM gcr.io/self-fluids/rocm:latest

COPY --from=devel /apps/self /apps/self
