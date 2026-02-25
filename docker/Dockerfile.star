FROM continuumio/miniconda3

LABEL description="STAR environment for RNA-seq alignment (STAR v2.7.11b)"
LABEL maintainer="ASE-network-clustering"

WORKDIR /app

# Copy environment file
COPY environments/star_env.yml .

# Create conda environment
RUN conda env create -f star_env.yml && \
    conda clean -afy

# Set environment as default
ENV PATH=/opt/conda/envs/star_env/bin:$PATH
SHELL ["conda", "run", "-n", "star_env", "/bin/bash", "-c"]

# Verify installation
RUN STAR --version

CMD ["/bin/bash"]
