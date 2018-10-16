# delphys

============================

## How to create or update BuildFile.xml
```
cd analysis/test
./makeBuildFile --d=/path/to/Delphes
```
=======

## How to git clone
```
scram p -n delphys CMSSW CMSSW_9_4_4
cd delphys/src
cmsenv
git clone git@github.com:"your_user_name"/nano.git #"your_user_name" should be replaced by real user name
scram b -j 20
```
## Download and untar delphes
```
cmsenv #in the delPhys directory
cd $CMSSW_BASE/src/nano/analysis/test/
mkdir hh
cd hh
wget http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.4.1.tar.gz
tar -xzvf Delphes-3.4.1.tar.gz
cd Delphes-3.4.1
sed -i s:c++0x:c++17: Makefile
make -j 4
export LD_LIBRARY_PATH=$CMSSW_BASE/src/nano/analysis/test/hh/Delphes-3.4.1:$LD_LIBRARY_PATH
```
## Add delphes library path into the .bashrc (or .bash_profile)
```
cat >> ~/.bashrc <<EOL
#Delphes Library Path
export LD_LIBRARY_PATH=$CMSSW_BASE/src/nano/analysis/test/hh/Delphes-3.4.1:\$LD_LIBRARY_PATH
EOL
```

## For doublehiggs analyser
```
please read delphys/analysis/test/hh/recipe_for_doublehiggs.sh to set external environments
```