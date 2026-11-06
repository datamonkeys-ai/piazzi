# Piazzi in a container: builds the C++ core and the CLI (tests do not run in the image build — use ctest
# locally) and serves the web version (WebAssembly prebuilt in web/).
#   docker build -t piazzi .
#   docker run --rm -p 8777:8777 piazzi            # web version at http://localhost:8777
#   docker run --rm -v "$PWD/imgs:/data" piazzi piazzi --help
# debian:bookworm-slim (pinned by digest; update with `docker buildx imagetools inspect debian:bookworm-slim`)
FROM debian@sha256:99c0b25c46d0214d5a2091092e7908907dcaecd6a962cde0ac700acee222b6f3 AS build
RUN apt-get update && apt-get install -y --no-install-recommends cmake ninja-build g++ && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY cpp/ cpp/
RUN cmake -S cpp -B cpp/build -G Ninja -DASTRO_GUI=OFF -DCMAKE_BUILD_TYPE=Release \
 && cmake --build cpp/build

# python:3.13-slim (pinned by digest; update with `docker buildx imagetools inspect python:3.13-slim`)
FROM python@sha256:2d44738ebf4ce5f63b9d9ea871e9c7cfdc0426f7b0c97cac915fd88ded1d86bb AS runtime
RUN useradd --create-home --uid 10001 astro
COPY --from=build /src/cpp/build/piazzi /usr/local/bin/piazzi
COPY web/ /app/web/
USER astro
WORKDIR /app
ENV PORT=8777 HOST=0.0.0.0 LOG_LEVEL=INFO
EXPOSE 8777
HEALTHCHECK --interval=30s --timeout=3s CMD python3 -c "import os,urllib.request; urllib.request.urlopen('http://127.0.0.1:%s/healthz' % os.environ['PORT'], timeout=2)"
CMD ["python3", "web/serve.py"]
