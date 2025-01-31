FROM ubuntu:22.04

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive TZ=Europe/London apt-get -y install tzdata
RUN apt-get install -y gcc g++ curl libssl-dev build-essential pkg-config git python3 python3-dev python3-pip libhdf5-dev
RUN pip3 install rdkit pybind11 numpy scipy pkgconfig

WORKDIR /workspace

CMD ["echo", "testing"]