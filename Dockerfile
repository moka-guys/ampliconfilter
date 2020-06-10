FROM python:3.7
RUN mkdir /code /sandbox /resources
WORKDIR /code
ADD . /code
RUN pip3 install pysam==0.15.4 cython==0.29.19 

ENTRYPOINT ["/code/ampliconfilter/ampliconFilter.py"]

