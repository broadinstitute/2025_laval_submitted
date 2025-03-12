#!/bin/bash

curl -O https://zenodo.org/records/15000710/files/raw_FACS.tar.gz

mkdir -p raw_inputs

tar -xvf raw_FACS.tar.gz -C raw_inputs

rm raw_FACS.tar.gz