# contextbook for the TEAPOT pipeline
FROM eclipse-temurin:17-jre-jammy

RUN apt-get update \
 && apt-get install -y --no-install-recommends libgomp1 procps \
 && rm -rf /var/lib/apt/lists/*

ARG VERSION
ENV JAVA_OPTS=""

COPY "target/encyclopedia-${VERSION}-executable.jar" /opt/encyclopedia/encyclopedia.jar

WORKDIR /data

ENTRYPOINT ["/bin/sh", "-c", \
  "exec java $JAVA_OPTS -Djava.awt.headless=true -cp /opt/encyclopedia/encyclopedia.jar \"$@\"", "--"]

CMD ["edu.washington.gs.maccoss.encyclopedia.Encyclopedia", "-h"]
