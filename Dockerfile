
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
