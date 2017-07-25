#!/bin/bash

# Copy over autogrow source code
if [ -f autogrow ]; then
    rm -r autogrow
fi
cp -r ../autogrow ./

# Build the autogrow docker image
docker build -t autogrow .

# Remove copied autogrow directory
rm -r autogrow