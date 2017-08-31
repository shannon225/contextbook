#!/bin/sh

if [ "$#" -ne 2 ]; then
	echo "You must specify a current and a target version (e.g. \"increment_version.sh 0.4.7 0.4.8\")";
	exit 1;
fi

echo "Updating from [${1}] to [${2}]";

mvn versions:set -DnewVersion="${1}"
hg commit -m "Update to version for release thesaurus-${1}."
hg tag "thesaurus-${1}" # optional
mvn clean package;
mvn versions:set -DnewVersion="${2}-SNAPSHOT"
hg commit -m "Update to next SNAPSHOT version."

echo "Finished updating from [${1}] to [${2}]";
