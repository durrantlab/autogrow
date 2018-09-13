#!/usr/bin/env python

import os
import sys
import shutil

# Check and make sure user input is alright.
sys_argv_lc = [a.lower() for a in sys.argv]
for prohibited_args in ["-vina_executable", "-openbabel_bin_directory", "-mgltools_directory", "-output_dir"]:
    if prohibited_args in sys_argv_lc:
        print "ERROR: Do not specify " + prohibited_args + " when using docker image!"
        sys.exit(0)

# Make a temporary directory to store files.
if os.path.exists(".tmp"):
    shutil.rmtree(".tmp")
os.mkdir(".tmp")

# Copy the receptor and source-compound files
recep_filename = sys.argv[sys.argv.index("-filename_of_receptor") + 1]
directory_of_source_compounds = sys.argv[sys.argv.index("-directory_of_source_compounds") + 1]
shutil.copyfile(recep_filename, "./.tmp/receptor.pdb")
shutil.copytree(directory_of_source_compounds, "./.tmp/ligands")
full_tmp_dir = os.path.abspath("./.tmp/")

# What about the fragments?
frag_dir_index = sys.argv.index("-directory_of_fragments") + 1
if sys.argv[frag_dir_index] in ["MW_150", "MW_200", "MW_250"]:
    sys.argv[frag_dir_index] = "/autogrow/fragments/" + sys.argv[frag_dir_index] + "/"
else:
    shutil.copytree(sys.argv[frag_dir_index], "./.tmp/frags")
    sys.argv[frag_dir_index] = "/autogrow_work_dir/frags/"

# Change additional parameters
def change_parameter(param_name, new_val):
    index = sys.argv.index("-" + param_name) + 1
    sys.argv[index] = new_val

change_parameter("filename_of_receptor", "/autogrow_work_dir/receptor.pdb")
change_parameter("directory_of_source_compounds", "/autogrow_work_dir/ligands/")

# Specify output directory
sys.argv.append("-output_dir")
sys.argv.append("/autogrow_work_dir/autogrow_output/")

# Make sure parameters with spaces are quoted
for i, a in enumerate(sys.argv):
    if " " in a:
        sys.argv[i] = "'" + '"' + a + '"' + "'"

os.system("""
sudo docker run \
    --rm \
    -it `# interactive mode` \
    -v """ + full_tmp_dir + """:/autogrow_work_dir/ \
    --name autogrow `# so a random name is not assigned` \
    autogrow `# the name` \
    """ + " ".join(sys.argv[1:]) + """ `# all the parameters`
""".strip()
)

# Move output to cwd
shutil.move("./.tmp/autogrow_output.zip", "./autogrow_output.zip")

# Delete tmp files
shutil.rmtree("./.tmp/")

# -d `# detached` \