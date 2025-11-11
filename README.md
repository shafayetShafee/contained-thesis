# Reproducibility Materials

[![R version](https://img.shields.io/badge/R-4.3.1-blue)](https://www.r-project.org/)
[![Reproducible Environment: renv](https://img.shields.io/badge/reproducible%20environment-renv-lightgreen.svg)](https://rstudio.github.io/renv/)
[![Formatted with air v0.7.1](https://img.shields.io/badge/formatted%20with-air%200.7.1-purple)](https://github.com/posit-dev/air)
[![License: MIT](https://img.shields.io/badge/License-MIT-orange)](https://opensource.org/licenses/MIT)
[![Docker Image](https://img.shields.io/badge/docker-ghcr.io/shafayetshafee/contained--thesis-blue?logo=docker)](https://ghcr.io/shafayetshafee/contained-thesis)
[![Docker](https://github.com/shafayetShafee/contained-thesis/actions/workflows/docker-publish.yml/badge.svg)](https://github.com/shafayetShafee/contained-thesis/actions/workflows/docker-publish.yml)

This repository contains the code materials and docker image to reproduce the simulation findings of the manuscript, 

> On the estimation of the median odds ratio for measuring contextual eﬀects in multilevel binary data from complex survey designs


## Directories & files descriptions

```yaml
├── Dockerfile              # Docker image build configuration for reproducible runs
├── docker-compose.yml      # Compose setup for running the project
├── docker-compose-dev.yml  # Development / interactive workflow setup
├── create_rprofile.sh      # Script to generate .Rprofile inside container

├── R/                      # Core R scripts containing simulation functions and runner.
├── plots/                  # Generated plots and visualization outputs
├── sim-results/            # Stored simulation output results
├── log/                    # Log text that will store the iteration-wise simulation status.

├── renv/                   # Project-local R package library (managed by `renv`)
└── renv.lock               # Snapshotted R package versions for reproducibility

├── LICENSE.md              # Project license (MIT)
├── README.md               # Project documentation
└── mor-interval.Rproj      # RStudio project file
```

## Steps to reproduce the simulation results

There are two ways to reproduce the analysis:

1. Using [`renv`](https://rstudio.github.io/renv/articles/renv.html) [local installation]

2. Using [Docker Image](https://www.docker.com/) [recommended]


### Method 1 — Using renv (Local Setup)

**Prerequisites:**

- [Git](https://git-scm.com/)
- [R](https://www.r-project.org/) [version 4.3.1 recommended]
- [RStudio](https://posit.co/downloads/)

**Steps**

1. Clone the repository
   ```bash
   https://github.com/shafayetShafee/contained-thesis.git
   ```

2. Enter the project directory
   ```bash
   cd contained-thesis
   ```
  
3. Open the project in RStudio by clicking the file `mor-interval.Rproj`.

4. When the project loads, renv will automatically activate and install all
required packages. This may take a few minutes on first run.

5. After installation completes, you can re-run the analysis scripts located in
the directory `R/`. Specifically The `sims_*` files contain the code to reproduce the results.

**Note on Reproducibility with `renv`:**

This project was originally developed using R version `4.3.1`. If your system uses 
a different R version, some packages may fail to install or behave differently. Therefore, 
this method may not guarantee fully reproducible results. If exact reproducibility 
is required, please use the Docker method below.


### Method 2 — Using Docker (Recommended)

**Prerequisite:** Docker installed on your system.

**Steps:**

1. Run the container using the pre-built image from GitHub.  
   (This image includes the correct R version, RStudio Server, and all R package
   dependencies exactly as used in the analysis.)
   ```bash
   docker run -d -p 8787:8787 ghcr.io/shafayetshafee/contained-thesis:2.0.0
   ```

2. Open your web browser and go to: `http://localhost:8787`

3. You will be connected to an RStudio Server session with the project
   pre-loaded and all dependencies installed. You can directly run the simulation
   scripts inside the RStudio Server.