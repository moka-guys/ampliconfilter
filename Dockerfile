FROM python:buster
RUN apt-get update && apt-get install -y samtools
RUN pip3 install cython==0.29.24
RUN pip3 install pysam==0.15.4
RUN mkdir /code /sandbox /resources
WORKDIR /code
ADD . /code

# ENTRYPOINT ["/usr/local/bin/python"]
CMD ["./ampliconFilter.py"]

