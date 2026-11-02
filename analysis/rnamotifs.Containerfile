FROM docker.io/rocker/r-ver:4.5.2

RUN apt-get update \
    && apt-get install -y --no-install-recommends git python3 cmake g++ make gfortran libpng-dev libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*
RUN git clone https://github.com/matteocereda/RNAmotifs.git /opt/RNAmotifs \
    && git -C /opt/RNAmotifs checkout 94d46a76350ed8c7d78cbce111703e2ccbab4fe3 \
    && cmake -S /opt/RNAmotifs -B /opt/RNAmotifs/build \
    && cmake --build /opt/RNAmotifs/build --parallel 2
RUN Rscript --vanilla -e 'install.packages(c("lattice", "latticeExtra", "bootstrap", "ggplot2", "reshape2"), repos="https://cloud.r-project.org"); stopifnot(all(vapply(c("lattice", "latticeExtra", "bootstrap", "ggplot2", "reshape2"), requireNamespace, logical(1), quietly=TRUE)))'

ENV R_PROFILE_USER=/dev/null
WORKDIR /opt/RNAmotifs
ENTRYPOINT ["python3", "/opt/RNAmotifs/rnamotifs"]
