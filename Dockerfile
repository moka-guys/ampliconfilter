FROM python:3.7
RUN apt-get update && apt-get install -y samtools
RUN pip3 install pysam==0.15.4 cython==0.29.19 
RUN mkdir /code /sandbox /resources
WORKDIR /code
ADD ampliconfilter/ /code

ENTRYPOINT [""]

