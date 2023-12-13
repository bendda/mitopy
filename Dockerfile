FROM python:3.11-buster as poetry-builder

RUN pip install poetry==1.4.2

ENV POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1 \
    POETRY_CACHE_DIR=/tmp/poetry_cache

WORKDIR /app

COPY pyproject.toml poetry.lock ./
COPY . .

RUN --mount=type=cache,target=$POETRY_CACHE_DIR poetry install --no-root
RUN poetry run pip install .


FROM debian:bullseye-20211220-slim AS gatk4-build
ARG GATK_VERSION=4.4.0.0
ENV BUILD_PACKAGES unzip gcc
RUN apt-get update && apt-get install -y --no-install-recommends ${BUILD_PACKAGES}
WORKDIR /home
ADD https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip gatk.zip 
RUN unzip gatk.zip && \
    rm -rf gatk*/gatkdoc/ gatk*/gatk*spark.jar


FROM python:3.11-slim AS mitopy

ARG GATK_VERSION=4.4.0.0
ARG BWAMEM2_VERSION=2.2.1
ARG HAPLOCHECK_VERSION=1.3.3
ARG MOSDEPTH_VERSION=0.3.5
ARG HAPLOGREP3_VERSION=3.2.1

RUN mkdir /usr/local/bin/gatk4 /usr/local/bin/bwamem2

# Install OpenJDK-17
RUN apt-get update && \
    apt-get install -y openjdk-17-jre unzip lbzip2 && \
    apt-get clean;


ENV VIRTUAL_ENV=/app/.venv \
    PATH="/app/.venv/bin:$PATH"

COPY --from=poetry-builder ${VIRTUAL_ENV} ${VIRTUAL_ENV}

# haplocheck
ADD https://github.com/genepi/haplocheck/releases/download/v${HAPLOCHECK_VERSION}/haplocheck.zip haplocheck.zip
RUN unzip haplocheck.zip -d /usr/local/bin/haplocheck

# haplogrep3
ADD https://github.com/genepi/haplogrep3/releases/download/v${HAPLOGREP3_VERSION}/haplogrep3-${HAPLOGREP3_VERSION}-linux.zip haplogrep3.zip
RUN unzip haplogrep3.zip -d /usr/local/bin/haplogrep3

# mosdepth
ADD https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth /usr/local/bin/mosdepth
RUN chmod u+x /usr/local/bin/mosdepth

# bwamem2
ADD https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWAMEM2_VERSION}/bwa-mem2-${BWAMEM2_VERSION}_x64-linux.tar.bz2 bwa-mem2.tar.bz2
RUN tar -xf bwa-mem2.tar.bz2 -C /usr/local/bin/bwamem2 --strip-components=1 --no-same-owner

# snpEff
ADD https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip snpeff.zip
RUN unzip snpeff.zip -d /usr/local/bin 
RUN rm -rf /usr/local/bin/snpEff/examples/ /usr/local/bin/snpEff/galaxy/

RUN rm bwa-mem2.tar.bz2 haplocheck.zip snpeff.zip haplogrep3.zip

COPY --from=gatk4-build /home/gatk-${GATK_VERSION}/ /usr/local/bin/gatk4
COPY --from=gatk4-build /usr/lib/x86_64-linux-gnu/libgomp.so.1 /usr/lib/x86_64-linux-gnu
ENV PATH="/usr/local/bin/haplogrep3:/usr/local/bin/snpEff/exec/:/usr/local/bin/mosdepth:/usr/local/bin/haplocheck:/usr/local/bin/bwamem2:/usr/local/bin/gatk4:${PATH}"
CMD [ "bash" ]



