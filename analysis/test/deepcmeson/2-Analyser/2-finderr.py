import glob
import os,sys
import ROOT

''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''
'''                                    '''
'''  python 2-finderr.py [data_name]   '''
'''                                    '''
''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''

def main():
    # Set Input Folder Name 
    data_name = 'test'
    if len(sys.argv) == 2:
        data_name = sys.argv[1:][0]
    elif len(sys.argv) > 2: 
        print "Wrong Input"
        sys.exit()
    
    path_name = '/xrootd/store/user/yyoun/DeepCMeson/analyser_{}/*.root'.format(data_name)
    all_paths = glob.glob(path_name)
    num_total = len(all_paths)
    print path_name
    
    # Check If File is Zombie
    num_zombies = 0
    for path in all_paths:
        root_file = ROOT.TFile(path)
        if root_file.IsZombie():
            num_zombies += 1
            root_file.Close()
            del root_file

            print(path)
            xrd_path = os.path.join("/xrd/", path.lstrip("/xrootd/"))
            command = "xrdfs cms-xrdr.private.lo:2094 rm {}".format(xrd_path)
            os.system(command)

    print("Zombies / Total = {} / {} ({:.1f})".format(num_zombies, num_total, float(num_zombies) / num_total))

if __name__ == '__main__':
    main()
