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


docs="common logging tracing tasking hmemo kregistry lama"


#create all docs
for doc in $docs; do
    # create output folder
    mkdir output/$doc
    
    # copy all json files to the main folder
    find $1/$doc/doc/ \( -name "*.fjson" -o -name "*.json" \) -exec cp {} output/$doc \;
    
    # copy all additional files
    if [ -e $1/$doc/doc/_images ]; then
        cp -R $1/$doc/doc/_images output/$doc
    fi
    if [ -e $1/$doc/doc/_downloads ]; then
        cp -R $1/$doc/doc/_downloads output/$doc
    fi
    
    # rename main index file to fit the guidelines of the rest plugin
    mv output/$doc/index.fjson output/$doc/Index.fjson
    
    # fix download urls
    find output/$doc/* \( -name "*.fjson" -o -name "*.json" \) -exec sed -i 's/href=\\\"[\.\/]*_downloads/href=\\\"http:\/\/www.libama.org\/fileadmin\/LAMA\/docs\/$doc\/_downloads/g' {} \;

    # flatten all internal references to the top level (fixes problems with the rest plugin    
    find output/$doc/* \( -name "*.fjson" -o -name "*.json" \) -exec sed -i 's/class=\\\"reference internal\\\" href=\\\"\([/.]*\)[a-Z0-9/]*\/\([a-Z0-9+]\+\)\/\?\\\"/class=\\\"reference internal\\\" href=\\\"\1\2\\\"/g' {} \;

    # replace all inter-sphinx links with the corresponding links on the website
    for doc2 in $docs; do
        find output/$doc/* \( -name "*.fjson" -o -name "*.json" \) -exec sed -i 's/href=\\\"[\.\/0-9a-Z_-]*\/share\/doc\/scai-'"${doc2}"'-[0-9\.]*\/#main-page/href=\\\"http:\/\/www.libama.org\/documentation\/'"${doc2}"'.html/g' {} \;
    done
done