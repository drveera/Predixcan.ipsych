#!/bin/sh

chmod +x scripts/pipeline.sh
echo "alias predixcan='$PWD/scripts/pipeline.sh'" >> ~/.bashrc
source ~/.bashrc
