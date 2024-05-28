### Some instructions
# docker pull rocker/geospatial # get image
# touch Dockerfile # create docker file
# docker run -it rocker/geospatial sh # test things interactively in the docker
# docker build -t simbinomialsdm . # build docker image with app (use lower case)
# docker run -p 8180:8180 simbinomialsdm # run app from docker
# docker images # list images
# docker image # list processes
# docker image rm -f 812a84f022b2 # remove image by ID or name


# Base R Shiny image
FROM rocker/geospatial

# Make a directory in the container
RUN mkdir /home/shiny-app

# Install R dependencies
RUN R -e "options(repos = c( \
             INLA = 'https://inla.r-inla-download.org/R/testing', \
             CRAN = 'https://cloud.r-project.org' \
          )); \
          install.packages(c('shinyjs', 'sn', 'INLA', 'inlabru')); \
          remotes::install_github('daattali/shinycssloaders')"

# Copy the Shiny app code
COPY app.R /home/shiny-app/app.R

# Expose the application port
EXPOSE 8180

# Run the R Shiny app
CMD Rscript /home/shiny-app/app.R
