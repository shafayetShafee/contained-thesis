## How to Reproduce the Results of this Thesis

### Prerequisites

1.  [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
2.  [Docker](https://docs.docker.com/engine/install/)

### Steps

1.  Clone this repository,

    ``` bash
    git clone --depth 1 https://github.com/shafayetShafee/contained-thesis.git
    ```

    (This will download only the lasted version, without all the git history, which can get somehwat bloated.)

2.  Then change directory to the repository,

    ``` bash
    cd contained-thesis
    ```

3.  Then from the `contained-thesis` directory, run the following,

    ``` bash
    docker run --rm -p 8787:8787 -d  -e DISABLE_AUTH=true -v $(pwd):/home/rstudio/thesis -v /home/rstudio/thesis/renv kshafayet/contained-thesis
    ```

    Then open your web browser to [localhost:8787](http://127.0.0.1:8787/) and you'll be welcomed by an RStudio session with a thesis folder with all that you need.
    
4. Then open the `R` folder within the `thesis` folder from the `Files` tab in the opened Rstudio server in your browser.

5. The `sims_*` files contain the code to reproduce the results.
