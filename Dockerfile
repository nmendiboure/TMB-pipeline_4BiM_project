FROM debian:latest

RUN apt update
RUN apt -y upgrade 
RUN apt install -y python3-dev python3-pip perl build-essential libssl-dev libffi-dev

RUN pip3 install --upgrade pip
RUN pip3 install numpy
RUN pip3 install pandas

ADD ./WorkingDirectoryDocker .

ENTRYPOINT ["python3"]
CMD ["/src/main.py"]
