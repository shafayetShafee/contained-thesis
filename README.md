## How to Reproduce the Results of this Thesis

### Prerequisites

1. git 
2. Docker

### Steps

1. Clone this repository,

``` bash
git clone --depth 1 https://github.com/shafayetShafee/contained-thesis.git
```

(This will download only the lasted version, without all the git history, which can get somehwat bloated.)

2. Then change directory to the repository,

``` bash
cd contained-thesis
```

3. Then from the `contained-thesis` directory, run the following,

``` bash
docker run --rm -p 8787:8787 -d  -e DISABLE_AUTH=true -v $(pwd):/home/rstudio/thesis -v /home/rstudio/thesis/renv  kshafayet/contained-thesis
```