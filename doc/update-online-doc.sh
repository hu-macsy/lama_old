#!/bin/bash

SF_USER=jirikraus

DOC_DIR=${HOME}/workspace/LAMA/doc

DOXYGEN_BIN=/home/lama/doxygen-1.7.4/bin/doxygen

cd ${DOC_DIR}

${DOXYGEN_BIN} LAMA-web.Doxyfile

mv ${HOME}/workspace/LAMA/doc/doxygen-web/html ${HOME}/workspace/LAMA/doc/doxygen-web/doc

rsync -avP -e ssh ${DOC_DIR}/doxygen-web/doc ${SF_USER},libama@web.sourceforge.net:htdocs

mv ${HOME}/workspace/LAMA/doc/doxygen-web/doc ${HOME}/workspace/LAMA/doc/doxygen-web/html
