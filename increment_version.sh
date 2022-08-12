#!/bin/sh

program="encyclopedia";
#program="xcordia";:

# think about auto version finder here:
#minor=`awk '/<version>/ {split($1, a, "."); split(a[3], b, "-"); print b[1]; exit 1}' pom.xml`; 
#mp1=$((minor+1));

if [ "$#" -ne 2 ]; then
	echo "You must specify a current and a target version (e.g. \"increment_version.sh 0.4.7 0.4.8\")";
	exit 1;
fi

echo "Updating ${program} from [${1}] to [${2}]";
echo `find . -name '*.java' -exec wc {} \; | awk '{sum=sum+$1} END {print sum}'` "total lines of Java code"

TAG=${program}-${1}

mvn versions:set -DnewVersion="${1}"
git commit -am "Update to version for release ${program}-${1}."
git tag "${TAG}" # optional
mvn clean package -DbuildJars;
mvn versions:set -DnewVersion="${2}-SNAPSHOT"
git commit -am "Update to next SNAPSHOT version."

echo "Finished updating from [${1}] to [${2}]";
echo "Remember to run \`git push origin --tags\` or \`git push origin ${TAG}\`!"
