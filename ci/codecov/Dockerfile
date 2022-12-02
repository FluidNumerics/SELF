FROM ubuntu:20.04 as builder


ENV DEBIAN_FRONTEND=noninteractive   \
    LANGUAGE=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

RUN apt-get -yqq update \
 && apt-get -yqq install --no-install-recommends \
        wget \
 && rm -rf /var/lib/apt/lists/*

RUN wget --no-check-certificate https://github.com/codecov/uploader/releases/download/v0.3.2/codecov-linux && \
    mv codecov-linux /usr/local/bin/codecov && \
    chmod +x /usr/local/bin/codecov

FROM ubuntu:20.04
COPY --from=builder /usr/local/bin/codecov /usr/local/bin/codecov

ENTRYPOINT ["/usr/local/bin/codecov"]

