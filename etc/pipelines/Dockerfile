FROM ubuntu:16.04

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    libncurses-dev \
    unzip \
    zlib1g-dev

COPY cpp-std-namespace.patch /usr/local/src
RUN git clone https://github.com/statgen/topmed_freeze3_calling.git /usr/local/src/topmed_freeze3_calling
RUN cd /usr/local/src/topmed_freeze3_calling && patch -p1 < ../cpp-std-namespace.patch
RUN make -C /usr/local/src/topmed_freeze3_calling/htslib install
RUN make -C /usr/local/src/topmed_freeze3_calling
RUN ln -s /usr/local/src/topmed_freeze3_calling/vt/vt /usr/local/bin/vt
RUN ln -s /usr/local/src/topmed_freeze3_calling/samtools/samtools /usr/local/bin/samtools
RUN ln -s /usr/local/src/topmed_freeze3_calling/gotcloud/src/bamUtil/bin/bam /usr/local/bin/bam

RUN mkdir -p /usr/local/share/ref
COPY ref/hs38DH.fa* /usr/local/share/ref/
#RUN /usr/local/src/topmed_freeze3_calling/samtools/misc/seq_cache_populate.pl -root /usr/local/share/ref/md5 /usr/local/share/ref/hs38DH.fa
