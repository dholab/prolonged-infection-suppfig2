FROM continuumio/miniconda3

SHELL ["/bin/bash", "-c"]

ADD suppfig2.yaml /tmp/suppfig2.yaml

RUN conda env create -f /tmp/suppfig2.yaml -n suppfig2 python=3.8 \
	&& conda clean --all -f --yes \
	&& conda init bash \
	&& echo "conda activate suppfig2" > ~/.bashrc

ENV PATH /opt/conda/envs/suppfig2/bin:$PATH

ENV CONDA_DEFAULT_ENV suppfig2

CMD ["/bin/bash"]
