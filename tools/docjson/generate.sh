#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Invalid number of parameters!"
    echo "usage: generate.sh <path>"
    echo "where <path>: Build folder which holds all built scai subprojects with json docs"
    exit 1
fi

# clear and create new output folder
rm -rf output
mkdir output

# copy all files from the lama build folder into the output folder
find $1/lama/doc/ \( -name "*.fjson" -o -name "*.json" \) -exec cp {} output \;
# copy all lama images
cp -R $1/lama/doc/_images output
cp -R $1/lama/doc/_downloads output

# copy other documentations linked by intersphinx
cp -R $1/logging/doc/index.fjson output/scaiLogging.fjson
cp -R $1/tracing/doc/index.fjson output/scaiTracing.fjson



# post processing
mv output/index.fjson output/Index.fjson

find output/* \( -name "*.fjson" -o -name "*.json" \) -exec sed -i 's/href=\\\"[\.\/]*_downloads/href=\\\"http:\/\/www.libama.eu\/fileadmin\/LAMA\/json\/_downloads/g' {} \; 
find output/* \( -name "*.fjson" -o -name "*.json" \) -exec sed -i 's/href=\\\"[\.\/0-9a-Z_-]*\/share\/doc\/scai-logging-[0-9\.]*\/#main-page/href=\\\"scaiLogging\//g' {} \;
find output/* \( -name "*.fjson" -o -name "*.json" \) -exec sed -i 's/href=\\\"[\.\/0-9a-Z_-]*\/share\/doc\/scai-tracing-[0-9\.]*\/#main-page/href=\\\"scaiTracing\//g' {} \;
