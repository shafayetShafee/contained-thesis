# get the base image, the rocker/verse has R, RStudio and pandoc
FROM rocker/rstudio:4.3.1

# Get and install system dependencies
WORKDIR /home/rstudio/thesis
RUN sudo apt-get update -qq \
 && apt-get -y --no-install-recommends install libfontconfig1-dev libfreetype6-dev cmake make libpng-dev libicu-dev pandoc imagemagick libmagick++-dev gsfonts libxml2-dev libssl-dev libcurl4-openssl-dev


# Get and install R packages to local library
COPY renv.lock renv.lock
COPY renv/activate.R renv/activate.R
RUN echo "setwd(\"/home/rstudio/thesis/\")\nsource(\"renv/activate.R\")" > ~/../home/rstudio/.Rprofile
RUN chown -R rstudio . \
 && sudo -u rstudio R -e 'renv::restore()'
