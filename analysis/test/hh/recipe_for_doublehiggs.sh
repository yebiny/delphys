#### Set CMSSW Environment and Get Github Delphys Repository
#scram p -n delPhys CMSSW CMSSW_9_4_4
#cd delPhys/src
#cmsenv

cd $CMSSW_BASE/src
git clone git@github.com:"your_user_name"/delphys.git #"your_user_name" should be replaced by real user name
#build analysers 
scram b -j 20
getFiles

#### Delphes  
cd $CMSSW_BASE/src/delphys/external/interface/
wget http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.4.1.tar.gz
tar -xzvf Delphes-3.4.1.tar.gz
rm Delphes-3.4.1.tar.gz
cd Delphes-3.4.1
sed -i s:c++0x:c++17: Makefile #set c++ version as 17th in Makefile
make -j 4
export LD_LIBRARY_PATH=$CMSSW_BASE/src/delphys/external/interface/Delphes-3.4.1:$LD_LIBRARY_PATH
#add delphes library path into the .bashrc (or .bash_profile)
cat >> ~/.bashrc <<EOL
#Delphes Library Path
export LD_LIBRARY_PATH=$CMSSW_BASE/src/delphys/external/interface/Delphes-3.4.1:\$LD_LIBRARY_PATH
EOL

#### Oxbridge
cd $CMSSW_BASE/src/delphys/
wget http://www.hep.phy.cam.ac.uk/~lester/dtm662/mt2/Releases/oxbridgekinetics.tar.gz
tar -xzvf oxbridgekinetics.tar.gz
rm oxbridgekinetics.tar.gz
mv oxbridgekinetics-1.2/ oxbridgekinetics
cd oxbridgekinetics
touch BuildFile.xml
#####################################################################
#make a buildfile for integrating oxbridgekinetics library into nano
#specify flags in BuildFile to ignore unnecessary warnings and errors
cat > BuildFile.xml << EOL
<use name="root" />
<flags LDFLAGS="-lMinuit2" />
<flags CXXFLAGS="-Wno-error=unused-but-set-variable" />
<flags CXXFLAGS="-Wno-error=unused-variable" />
<flags CXXFLAGS="-Wno-error=sign-compare" />
<flags CXXFLAGS="-Wno-error=maybe-uninitialized" />
<export>
  <lib name="1"/>
</export>
EOL 
#####################################################################
mv Mt2/ src
cd src/Mt2
sed -i s:"Mt2/:": * # modify the path of Mt2 header files.
cd ../../
scram b
#set right path in the nano/analysis/bin/BuildFile.xml
sed -i s:CMSSW_BASE:$CMSSW_BASE: $CMSSW_BASE/src/delphys/analysis/bin/BuildFile.xml
sed -i '/erase/d' $CMSSW_BASE/src/delphys/analysis/bin/BuildFile.xml
cd $CMSSW_BASE/src/delphys/analysis/bin/
scram b -j12
doubleHiggsAnalyser
