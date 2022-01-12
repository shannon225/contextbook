FROM openjdk:8

RUN apt-get update && \
    apt-get -y upgrade && \ 
    apt-get -y install libgomp1 && \
    apt-get clean

WORKDIR /code

ARG VERSION=1.12.31
ENV VERSION ${VERSION}

RUN wget -O "encyclopedia-$VERSION-executable.jar" "https://bitbucket.org/searleb/encyclopedia/downloads/encyclopedia-$VERSION-executable.jar"

WORKDIR /app