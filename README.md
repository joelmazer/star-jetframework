# star-jetframework
This is a framework for:
1) jet-finding, 
2) event rho calculation
3) QA 

Additonally includes an event pool manager for mixed events, a correlations task and some base classes


This framework requires the installation of FastJet + FastJet contrib.
FastJet can be found at: http://fastjet.fr/all-releases.html
-download FastJet 3.0 +

FastJet contrib can be found at:
http://fastjet.hepforge.org/contrib/

RCAS & PDSF
PDSF requires 32-bit build to fit into other software framwork

FASTJET installation example:
  tar zxvf fastjet-3.2.2.tar.gz
  cd fastjet-3.2.2/
  ./configure --prefix=$PWD/../fastjet-install CXXFLAGS="-m32 -fPIC -fno-inline" CFLAGS="-m32 -fPIC -fno-inline" LDFLAGS="-m32"

  make
  make check
  make install


FastJet contrib installation example:
  tar zxvf fjcontrib-1.026.tar.gz
  cd fjcontrib-1.026
  ./configure --fastjet-config=$PWD/../fastjet-3.2.2/fastjet-config --prefix=$PWD/../fastjet-install CXXFLAGS="-m32 -fPIC -fno-inline" CFLAGS="-m32 -fPIC -fno-inline" LDFLAGS="-m32"

  make
  make install
  make fragile-shared
  make fragile-shared-install


unset flags after:
  unset CXXFLAGS CFLAGS LDFLAGS
