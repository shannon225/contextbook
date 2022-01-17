FROM openjdk:11.0-jre-slim

RUN apt-get update && \
    apt-get -y upgrade && \ 
    apt-get -y install libgomp1 wget && \
    apt-get clean

WORKDIR /code

ARG VERSION
ENV VERSION ${VERSION}

RUN wget -O "encyclopedia-$VERSION-executable.jar" "https://bitbucket.org/searleb/encyclopedia/downloads/encyclopedia-$VERSION-executable.jar"

WORKDIR /app