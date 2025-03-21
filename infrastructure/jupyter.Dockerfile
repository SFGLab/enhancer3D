FROM quay.io/jupyter/base-notebook:python-3.12.7 as jupyter

USER root

# Install Kerberos
RUN conda install requests-kerberos -y

# Install git
RUN apt-get update && apt-get install -y git zip unzip
RUN git config --global user.email "jupyter@sandbox.kot.tools"
RUN git config --global user.name "Jupyter"

USER $NB_USER

# Install Sparkmagic
RUN pip install --upgrade pip
RUN pip install --upgrade --ignore-installed setuptools
RUN pip install sparkmagic ipywidgets

RUN mkdir /home/$NB_USER/.sparkmagic
COPY infrastructure/config/jupyter/config.json /home/$NB_USER/.sparkmagic/config.json

RUN jupyter-kernelspec install --user $(pip show sparkmagic | grep Location | cut -d" " -f2)/sparkmagic/kernels/sparkkernel
RUN jupyter-kernelspec install --user $(pip show sparkmagic | grep Location | cut -d" " -f2)/sparkmagic/kernels/pysparkkernel
RUN jupyter-kernelspec install --user $(pip show sparkmagic | grep Location | cut -d" " -f2)/sparkmagic/kernels/sparkrkernel
RUN jupyter server extension enable --py sparkmagic

USER root

RUN chown $NB_USER /home/$NB_USER/.sparkmagic/config.json

USER $NB_USER

CMD ["start-notebook.sh", "--NotebookApp.iopub_data_rate_limit=1000000000"]
