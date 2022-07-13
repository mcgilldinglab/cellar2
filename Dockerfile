FROM euxhen/bioenv
LABEL author="ehasanaj@cs.cmu.edu"

RUN mkdir /home/nonroot/cellar
ARG VER=unknown
COPY . /home/nonroot/cellar
RUN sudo chown nonroot:nonroot /home/nonroot/cellar/data/uploaded
RUN sudo chown nonroot:nonroot /home/nonroot/cellar/tmp

WORKDIR /home/nonroot/cellar
EXPOSE 8050
# ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "cellar", "gunicorn", "-t", "7200", "--threads", "8", "-b", "0.0.0.0:8050", "main:server"]
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "cellar", "python", "main.py"]