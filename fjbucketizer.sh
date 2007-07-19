#!/bin/sh

DIR=`dirname $0`
. ${DIR}/Configuration
${PYTHON} ${FJBUCKETIZER}/fjbucketizer.py $*
