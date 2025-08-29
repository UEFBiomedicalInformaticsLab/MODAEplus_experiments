# Docker
## Build image

```
export IMAGE_NAME=MODAEplus
export PATH_TO_MODEL=../../MODAEplus
docker buildx build --build-context modae=$PATH_TO_MODEL --tag=$IMAGE_NAME
```

## Run container

```
export IMAGE_NAME=MODAEplus
export PATH_TO_SCRIPTS=../../MODAEplus_experiments/
export PATH_TO_DATA=../../MODAEplus_data/
docker run --rm \
  --mount type=bind,src=$PATH_TO_DATA,dst=/data \
  --mount type=bind,src=$PATH_TO_SCRIPTS,dst=/code \
  -it --net=host --env DISPLAY=$DISPLAY $IMAGE_NAME
```
