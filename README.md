# Requirements

Only a working ROOT installation is required.

# Build / Installation

To build and install the repository:

```bash
git clone https://baltig.infn.it/LUNA_DAQ/LFReader.git
cd LFReader && mkdir build && cd build/
cmake ..
make && sudo make install
```

The binaries will be installed in ```/opt/LFReader/```. You can add the ```export PATH=/opt/LFReader/${PATH}``` line to your ```bashrc``` file to start it from 
terminal. As an alternative you can create a symlink:

```bash
ln -s /opt/LFReader/LFReader /usr/bin/LFReader
```

## Usage
It is necessary to ID and channels of the boards that were used during the acquisition. For show help for the input arguments use:

```bash
LFReader --help
```
