FROM nvcr.io/nvidia/tensorflow:25.02-tf2-py3
#FROM tensorflow/tensorflow:2.19.0

RUN apt update && apt install -y spyder
COPY --from=modae ./ ./MODAE
RUN pip install --no-cache-dir ./MODAE
RUN mkdir /data
RUN mkdir /scripts

CMD spyder
