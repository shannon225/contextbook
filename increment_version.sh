#!/bin/sh

program="encyclopedia";
#program="xcordia";:

# think about auto version finder here:
#minor=`awk '/<version>/ {split($1, a, "."); split(a[3], b, "-"); print b[1]; exit 1}' pom.xml`; 
#mp1=$((minor+1));

#if [ "$#" -ne 2 ]; then
#	echo "You must specify a current and a target version (e.g. \"increment_version.sh 0.4.7 0.4.8\")";
#	exit 1;
#fi

branch=$(git symbolic-ref --short HEAD);
if [ $branch == "master" ]; then
	branch="develop";
fi

let "year=`date +%Y`-2020";
current=${year}.`date +%-m.%-d`;
let "next=${year}+1";
next=${next}.0.0;

echo "Updating ${program} $branch from [${current}] to [${next}]";
echo `find . -name '*.java' -exec wc {} \; | awk '{sum=sum+$1} END {print sum}'` "total lines of Java code"

TAG=${program}-${branch}-${current}

mvn versions:set -DnewVersion="${current}"
git commit -am "Update to version for release ${program}-${branch}-${current}."
git tag "${TAG}" # optional
mvn -Djavacpp.platform.custom -Djavacpp.platform.host \
	-Djavacpp.platform.macosx-x86_64 \
	-Djavacpp.platform.macosx-arm64 \
	-Djavacpp.platform.linux-x86_64 \
	-Djavacpp.platform.windows-x86_64 \
	clean package -DbuildJars;
mvn versions:set -DnewVersion="${next}-${branch}-SNAPSHOT"
git commit -am "Update to next SNAPSHOT version."

echo "Finished updating from [${current}] to [${next}]";
echo "Remember to run \`git push origin --tags\` or \`git push origin ${TAG}\`!"
