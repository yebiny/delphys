import glob
import os, sys
import ROOT

''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''
'''                                    '''
'''   python 3-mergy.py [data_name]    '''
'''                                    '''
''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''

def main():
    # Set Merge Input Folder 
    data_name = 'test'
    if len(sys.argv) == 2:
        data_name = sys.argv[1:][0]
    elif len(sys.argv) > 2 : 
        print "Wrong Input"
        sys.exit()
    
    path_name = '/xrootd/store/user/yyoun/DeepCMeson/analyser_{}/*.root'.format(data_name)
    all_paths = glob.glob(path_name)
    merge_dir = './result/{}'.format(data_name)
   
    # Make Merge Result Folder
    if not os.path.isdir(merge_dir):
        print "Make Folder: ", merge_dir
        os.mkdir(merge_dir)
    else: 
        print"Folder is already exist: ", merge_dir
        sys.exit()
   
    # Start Merge 500 pieces     
    names_target = ""
    for i, path in enumerate(all_paths):
        names_target += "%s "%(path)

        if (i+1)%500 == 0:
            command = 'hadd {}/deepc_{}.root {}'.format(merge_dir, i, names_target)
            os.system(command)
            # Reset merge target 
            names_target = ""

    command = 'hadd {}/deepc_{}.root {}'.format(merge_dir,i, names_target)
    os.system(command)

if __name__ == '__main__':
    main()
