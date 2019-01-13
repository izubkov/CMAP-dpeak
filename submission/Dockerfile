FROM izubkov/cmap:0

RUN mkdir -p /submission
WORKDIR /submission
COPY . /submission

RUN chmod +x /submission/exec_R.sh

ENTRYPOINT [ "/submission/exec_R.sh" ]
CMD []
