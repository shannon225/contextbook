#!/bin/sh

# think about auto version finder here:
#minor=`awk '/<version>/ {split($1, a, "."); split(a[3], b, "-"); print b[1]; exit 1}' pom.xml`; 
#mp1=$((minor+1));

if [ "$#" -ne 2 ]; then
	echo "You must specify a current and a target version (e.g. \"increment_version.sh 0.4.7 0.4.8\")";
	exit 1;
fi

echo "Updating from [${1}] to [${2}]";

mvn versions:set -DnewVersion="${1}"
hg commit -m "Update to version for release encyclopedia-${1}."
hg tag "encyclopedia-${1}" # optional
mvn clean package;
mvn versions:set -DnewVersion="${2}-SNAPSHOT"
hg commit -m "Update to next SNAPSHOT version."

echo "Finished updating from [${1}] to [${2}]";
