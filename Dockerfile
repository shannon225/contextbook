FROM openjdk:8-jre

RUN apt-get update && \
    apt-get -y upgrade && \ 
    apt-get -y install libgomp1 && \
    apt-get clean

WORKDIR /code

ARG VERSION
ENV VERSION ${VERSION}

# Copy the executable JAR from the build context to the container.
# Requires that a correctly-versioned JAR has been built locally.
# If you encounter an error here running `docker build` check that
# you've run `mvn package -PbuildEncyclopedia` and specified the
# correct `--build-arg VERSION`.
COPY "target/encyclopedia-$VERSION-executable.jar" .

WORKDIR /app

# Entry point for `docker run` specifies running the JAR in headless mode.
# We must use this combination of shell and exec so that variable expansion
# functions, parameters may be specified by CMD/run, and the process receives
# signals (e.g. SIGINT in case of ctrl-C or `docker stop`).
# See https://stackoverflow.com/a/42332740 for more information.
ENTRYPOINT ["/bin/sh", "-c", "exec java -jar /code/encyclopedia-$VERSION-executable.jar -Djava.awt.headless \"$0\" \"$@\""]

# CMD provides default parameters
CMD ["-h"]