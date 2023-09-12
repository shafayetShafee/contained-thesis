## TODOS

See the tasks in [#1](https://github.com/shafayetShafee/thesis/issues/1)

``` bash
docker run --rm -p 8787:8787 -d  -e DISABLE_AUTH=true -v $(pwd):/home/rstudio/thesis -v /home/rstudio/thesis/renv  kshafayet/contained-thesis
```