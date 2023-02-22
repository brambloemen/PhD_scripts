#!usr/bin/bash

# Set the filetype to be gzipped as the first argument
filetype=$1

# Find all files with the specified filetype and gzip them
find . -name "*.$filetype" -exec pigz -p9 {} \;