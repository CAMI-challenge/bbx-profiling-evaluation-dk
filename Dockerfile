FROM ubuntu:14.04
MAINTAINER Stefan Janssen, sjanssen@techfak.uni-bielefeld.de

RUN apt-get update -y
RUN apt-get install -y python2.7 python-numpy python-networkx python-scipy

ENV PREFIX /biobox
ENV OUTPUTDIR /output

ADD EMDUnifrac.py ${PREFIX}/bin/
ADD ProfilingMetrics.py ${PREFIX}/bin/

RUN chmod u+x ${PREFIX}/bin/EMDUnifrac.py ${PREFIX}/bin/ProfilingMetrics.py ${PREFIX}/bin/Utils.py

#ENTRYPOINT ["python", "/usr/local/bin/ProfilingMetrics.py", "-g", "/input_truth.tsv", "-r", "/input_reconstruction.tsv", "-o", "/output/output.txt", "-u", "/usr/local/bin/EMDUnifrac.py", "-e", "0", "-n", "y"]
