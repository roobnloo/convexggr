FROM rocker/verse:4.2

### install python
# copy from jupyter stacks
# roughly around the same time as rocker/verse:4.2
COPY --from=jupyter/base-notebook:x86_64-python-3.10 /opt/conda /opt/conda
ENV PATH /opt/conda/bin:$PATH

### install quarto
# versions: https://quarto.org/docs/download/_download.json
RUN curl --silent -L --fail \
        https://github.com/quarto-dev/quarto-cli/releases/download/v1.2.269/quarto-1.2.269-linux-amd64.deb > /tmp/quarto.deb && \
    apt update && \
    apt install -y --no-install-recommends /tmp/quarto.deb && \
    apt install -y --no-install-recommends tmux

USER rstudio

RUN pip install radian && \
    install2.r languageserver httpgd && \
    installGithub.r ManuelHentschel/vscDebugger


RUN install2.r abind RcppEigen sparsegl patchwork tictoc CVXR tables kableExtra

USER root