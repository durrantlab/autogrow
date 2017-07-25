# Pull base image.
FROM ubuntu:16.04

MAINTAINER Jacob Durrant <durrantj@pitt.edu>

# Install programs through apt-get
RUN apt-get -yqq update
RUN apt-get -fyqq install
RUN apt-get -yqq install wget
RUN apt-get -yqq update
RUN apt-get -yqq install openbabel
RUN apt-get -yqq install zip

# Install python
RUN apt-get -yqq install bzip2
RUN wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
RUN chmod +x Miniconda2-latest-Linux-x86_64.sh
RUN ./Miniconda2-latest-Linux-x86_64.sh -b
RUN rm ./Miniconda2-latest-Linux-x86_64.sh
RUN echo "alias python=/root/miniconda2/bin/python" >> /root/.bashrc

# Install python modules required to run AutoGrow
RUN /root/miniconda2/bin/conda install -y scipy

# Install mgltools
RUN wget http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/mgltools_x86_64Linux2_1.5.6.tar.gz
RUN tar xvfz *.tar.gz
RUN rm *.tar.gz
RUN cd /mgltool*/ && \
    ./install.sh
# RUN echo "alias python=/mgltools*/bin/pythonsh" >> /root/.bashrc

# Install AutoDock Vina
RUN wget http://vina.scripps.edu/download/autodock_vina_1_1_2_linux_x86.tgz
RUN tar xvfz *.tgz
RUN rm *.tgz

# Copy over autogrow files
RUN mkdir /autogrow/
ADD ./autogrow /autogrow

# Modify those files given current file system
RUN export mgldir=`ls -d /mgltools*/` && \
    export vina_bin=`ls /*vina*/bin/vina` && \
    cd /autogrow/ && \
    cat autogrow_3_1_2.py | sed "s|vars\['mgltools_directory'\] = \"\"|vars['mgltools_directory'] = \"${mgldir}\"|g" > tmp && \
    mv tmp autogrow_3_1_2.py && \
    cat autogrow_3_1_2.py | sed "s|vars\['openbabel_bin_directory'\] = \"\"|vars['openbabel_bin_directory'] = \"/usr/bin/\"|g" > tmp && \
    mv tmp autogrow_3_1_2.py && \
    cat autogrow_3_1_2.py | sed "s|vars\['vina_executable'\] = \"\"|vars['vina_executable'] = \"${vina_bin}\"|g" > tmp && \
    mv tmp autogrow_3_1_2.py && \
    cat autogrow_3_1_2.py | sed "s|vars\['babel_executable'\] = \"\"|vars['babel_executable'] = \"/usr/bin/babel\"|g" > tmp && \
    mv tmp autogrow_3_1_2.py && \
    cat autogrow_3_1_2.py | sed "s|vars\['obprop_executable'\] = \"\"|vars['obprop_executable'] = \"/usr/bin/obprop\"|g" > tmp && \
    mv tmp autogrow_3_1_2.py && \
    cat autogrow_3_1_2.py | sed "s|vars\['prepare_ligand4.py'\] = \"\"|vars['prepare_ligand4.py'] = \"${mgldir}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py\"|g" > tmp && \
    mv tmp autogrow_3_1_2.py && \
    cat autogrow_3_1_2.py | sed "s|vars\['prepare_receptor4.py'\] = \"\"|vars['prepare_receptor4.py'] = \"${mgldir}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py\"|g" > tmp && \
    mv tmp autogrow_3_1_2.py && \
    cat autogrow_3_1_2.py | sed "s|vars\['mgl_python'\] = \"\"|vars['mgl_python'] = \"${mgldir}/bin/pythonsh\"|g" > tmp && \
    mv tmp autogrow_3_1_2.py

# Copy over the autogrow run script
ADD ./run_autogrow_in_container.bash /autogrow/run_autogrow_in_container.bash
RUN chmod +x /autogrow/run_autogrow_in_container.bash

# Make directories for mounted host file systems
RUN mkdir /autogrow_work_dir

ENTRYPOINT ["bash", "/autogrow/run_autogrow_in_container.bash"]
