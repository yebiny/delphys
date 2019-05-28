#!/usr/bin/env python
import os,sys
import pprint

data_name = sys.argv[1]
pp = pprint.pprint
filename = './txtFiles/deepc_{}.txt'.format(data_name)
f = open(filename , "w")
a = []
pathlist=  [
            #"/xrootd/store/user/seyang/DeepCMeson/1-Delphes/"            
            #"/xrootd/store/user/seya-ng/TTvsQCD/1-Delphes/TTJets_aMC"     
            '/xrootd/store/user/yyoun/DeepCMeson/2-Analyser/'+data_name
]   
for path in pathlist:
    temp = os.listdir(path)
    temp = map(lambda p: os.path.join(path, p), temp) 
    a +=  temp
for i in a:
    newstr = str.replace(i, '/xrootd','root://cms-xrdr.private.lo:2094///xrd'
          # ' root://cms-xrdr.sdfarm.kr:1094///xrd'
            )
    f.write(newstr + "\n")

pp(newstr)
f.close()    
