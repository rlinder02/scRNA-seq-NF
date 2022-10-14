FROM redhat/ubi8:latest

ARG star_version=2.7.10a_alpha_220818

COPY conda.yml .

#Install OS packages
RUN yum -y install gcc gcc-c++ make

RUN curl -sSL https://repo.continuum.io/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local/ \
    && rm -rf /tmp/miniconda.sh

RUN conda env update -n root -f conda.yml && conda clean -a

#Install STAR
RUN wget --no-check-certificate https://github.com/alexdobin/STAR/archive/refs/tags/${star_version}.tar.gz && \
    tar -xzf ${star_version}.tar.gz -C /opt && \
    cd /opt/STAR-${star_version}/source && \
    make STAR CXXFLAGS_SIMD="-msse4.2" && \
    cd / && rm ${star_version}.tar.gz

ENV PATH="${PATH}:/opt/STAR-${star_version}/source" 






