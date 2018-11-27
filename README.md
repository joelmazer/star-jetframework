# star-jetframework
*this will be updated...*
This is a framework for:
* a FastJet wrapper, with various functionality 
* jet-finding
* underlying event rho calculation
* QA 
* creating reduced track ojbects (FemtoTracks)
* event mixing (for acceptance and detector effects via event pool manager)
* a base class use for inheritance
* studyng jet shapes
* calculating the event plane using multiple detectors + all included corrections
* jet-hadron corrrelations
* dummy analysis class to get started


Additionally contained are multiple other classes useful for the main correlation analysis

# Getting started
## Prequisites / also needed
inside your StRoot director you will need to copy the contents of this git repository into the directory
```
StMyAnalysisMaker/
```

For centrality, you will need the repository: 
```
https://github.com/joelmazer/Centrality
```
In particular, you will need the StRefMultCorr directory in the Centrality repository to be placed into your StRoot directory. This leads to the following StRoot structure:
```
StRoot/StMyAnalysisMaker/
StRoot/StRefMultCorr/
```

## Installing

### FASTJET install 
This framework requires the installation of FastJet + FastJet contrib.
FastJet can be found at: http://fastjet.fr/all-releases.html
* download FastJet 3.0 +

FastJet contrib can be found at:
```
http://fastjet.hepforge.org/contrib/
```

RCAS & PDSF both require 32-bit build to fit into other software framework.  To force this install on the 64-bit systems, please follow the below prescription to avoid headaches.

FASTJET installation example:
```
  tar zxvf fastjet-3.2.2.tar.gz
  cd fastjet-3.2.2/
  ./configure --prefix=$PWD/../fastjet-install CXXFLAGS="-m32 -fPIC -fno-inline" CFLAGS="-m32 -fPIC -fno-inline" LDFLAGS="-m32"

  make
  make check
  make install
```

FastJet contrib installation example:
```
  tar zxvf fjcontrib-1.026.tar.gz
  cd fjcontrib-1.026
  ./configure --fastjet-config=$PWD/../fastjet-3.2.2/fastjet-config --prefix=$PWD/../fastjet-install CXXFLAGS="-m32 -fPIC -fno-inline" CFLAGS="-m32 -fPIC -fno-inline" LDFLAGS="-m32"

  make
  make install
  make fragile-shared
  make fragile-shared-install
```

unset flags after:
```
  unset CXXFLAGS CFLAGS LDFLAGS
```

You will additionally need to set FastJet environmental variables in your startup file (mine is .bashrc, and the below is an example):
```
export FASTJET='/path/to/your/FastJet/fastjet-install'
export FASTJET_CONTRIB='/path/to/your/FastJet/fjcontrib-1.026'
```
and soft link the include paths to your working directory which contains your StRoot directory and where you run your readPicoDst.C file from:
```
ln -s /path/to/your/FastJet/include/siscone siscone
ln -s /path/to/your/FastJet/include/fastjet fastjet
```

## Build procedures
Once FastJet is properly installed


## Class descriptions
*Will be updated*

## Author
**Joel Mazer**

## Notes

