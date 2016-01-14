FROM ubuntu:14.04
MAINTAINER Stefan Janssen, sjanssen@techfak.uni-bielefeld.de

RUN apt-get update -y
RUN apt-get install -y python2.7 python-numpy python-networkx python-scipy wget ca-certificates xz-utils

ENV PREFIX /usr/local
ENV OUTPUTDIR /output

# Locations for biobox file validator
ENV VALIDATOR /bbx/validator/
ENV BASE_URL https://s3-us-west-1.amazonaws.com/bioboxes-tools/validate-biobox-file
ENV VERSION  0.x.y
RUN mkdir -p ${VALIDATOR}

# download the validate-biobox-file binary and extract it to the directory $VALIDATOR
RUN wget \
      --quiet \
      --output-document -\
      ${BASE_URL}/${VERSION}/validate-biobox-file.tar.xz \
    | tar xJf - \
      --directory ${VALIDATOR} \
      --strip-components=1

ENV PATH ${PATH}:${VALIDATOR}

#~ # download the assembler schema
#~ RUN wget \
    #~ --output-document /schema.yaml \
    #~ https://raw.githubusercontent.com/bioboxes/rfc/master/container/short-read-assembler/input_schema.yaml

ENV CONVERT https://github.com/bronze1man/yaml2json/raw/master/builds/linux_386/yaml2json
# download yaml2json and make it executable
RUN cd /usr/local/bin && wget --quiet ${CONVERT} && chmod 700 yaml2json

ENV JQ http://stedolan.github.io/jq/download/linux64/jq
# download jq and make it executable
RUN cd /usr/local/bin && wget --quiet ${JQ} && chmod 700 jq

ADD EMDUnifrac.py ${PREFIX}/bin/
ADD ProfilingMetrics.py ${PREFIX}/bin/

RUN chmod u+x ${PREFIX}/bin/EMDUnifrac.py ${PREFIX}/bin/ProfilingMetrics.py

# Add evaluate script to the directory ${PREFIX}/bin inside the container.
# ${PREFIX} is appended to the $PATH variable what means that every script
# in that directory will be executed in the shell  without providing the path.
ADD evaluate ${PREFIX}/bin/
ENV PATH ${PATH}:${PREFIX}/bin
RUN chmod -R u+x ${PREFIX}/bin
ADD Taskfile /

ENTRYPOINT ["evaluate"]
