FROM openanalytics/r-base


# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libssl1.1

USER root

# basic shiny functionality
RUN R -e "install.packages(c('shiny', 'rmarkdown'), repos='https://cloud.r-project.org/')"

# install dependencies of the multiome app
RUN R -e "install.packages(c('pacman', 'visNetwork', 'ggiraph'))"
RUN sudo apt install -y libcurl4-openssl-dev libxml2-dev
RUN R -e "install.packages(c('xml2', 'rvest'))"
RUN R -e "install.packages('tidyverse')"
RUN R -e "install.packages('cowplot')"
RUN R -e "install.packages('igraph')"

# copy the app to the image
RUN mkdir /root/Explorer
COPY Explorer /root/Explorer

COPY deploy/Rprofile.site /usr/lib/R/etc/

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/root/Explorer')"]
