#!/bin/bash

# Change into autogrow directory
cd /autogrow/

# Run autogrow
/root/miniconda2/bin/python ./autogrow_3_1_2.py $@

# Zip up results and move to working directory.
cd /autogrow_work_dir/
zip -r autogrow_output.zip autogrow_output/

# For interactive
bash