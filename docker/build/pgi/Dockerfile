FROM gcr.io/self-fluids/pgi as devel 
ARG BUILD_FLAGS="--enable-diagnostics --enable-tecplot"

ENV FCFLAGS="-O3" \
    cuda_toolkit_version="cuda10.1"

RUN mkdir -p /build

COPY . /build/self-fluids

RUN cd /build/self-fluids && \
    CPP=cpp \
    /build/self-fluids/configure --with-metis=${LIB_METIS}\
                                 ${BUILD_FLAGS} \
                                 --prefix=/usr/local/self-fluids || cat config.log && \
    make && make install


RUN mkdir -p /examples && \
    cp -r /build/self-fluids/examples/* /examples/

FROM gcr.io/self-fluids/pgi

 
COPY --from=devel /usr/local/self-fluids /usr/local/self-fluids

COPY --from=devel /examples /examples

