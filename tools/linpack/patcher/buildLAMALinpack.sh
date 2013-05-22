#!/bin/bash

PATCH_FILE=$LAMA_HOME/tools/linpack/patcher/lama.patch

if [ $# -eq 1 ]; then
	if [ -f $1 ]; then
		PATCH_FILE=$1
	fi
fi

if [ -f $PATCH_FILE ]; then
    echo "Using patch file $PATCH_FILE"
else
    echo "Could not found the patch file at $PATCH_FILE, please specify as first parameter"
fi

wget http://www.netlib.org/benchmark/hpl/hpl-2.0.tar.gz

tar -xzf hpl-2.0.tar.gz

cd hpl-2.0

patch -p1 < $PATCH_FILE