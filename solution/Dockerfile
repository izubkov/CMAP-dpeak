# Dockerized version of DECONV > DE transformation code, provided by Harvard.
FROM cmap/sig_2to4_tool

RUN apt update && apt install -y apt-transport-https software-properties-common
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu trusty-cran35/" && apt update && apt install --force-yes -y r-base

# install cmapR package

RUN Rscript -e "install.packages(c(\"BiocManager\", \"data.table\"))"
RUN Rscript -e "BiocManager::install(\"dplyr\")"
RUN Rscript -e "BiocManager::install(\"pROC\")"
RUN Rscript -e "BiocManager::install(\"rhdf5\")"
COPY ./cmapR_1.0.1.tar.gz /score/cmapR_1.0.1.tar.gz
RUN Rscript -e "install.packages(\"/score/cmapR_1.0.1.tar.gz\", type=\"source\", repos=NULL)"

# install additional solution packages

RUN Rscript -e "install.packages(\"magrittr\")"
RUN Rscript -e "install.packages(\"stringr\")"
RUN Rscript -e "install.packages(\"getopt\")"
RUN Rscript -e "install.packages(\"readr\")"
RUN Rscript -e "install.packages(\"purrr\")"

COPY . /solution

ENTRYPOINT [ "/solution/exec_R.sh" ]
