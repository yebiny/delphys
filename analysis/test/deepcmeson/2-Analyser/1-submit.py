import os, sys, getopt

data_name = 'test'
txt_name = '../1-Generation/txtFiles/deepcmeson_pwg_test.txt'
if len(sys.argv) == 2 :
    data_name = sys.argv[1:][0]
elif len(sys.argv) == 3 :
    data_name = sys.argv[1:][0]
    txt_name = sys.argv[1:][1]

_JDS_FMT = """executable = /cms/ldap_home/yyoun/delPhys/bin/slc6_amd64_gcc630/analyseDeepCMeson
universe   = vanilla

arguments = {in_path} {out_path} 

log = condor.log

requirements = ( HasSingularity == true )
accounting_group=group_cms
+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el6:latest"
+SingularityBind = "/cvmfs, /cms, /share"

getenv     = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = ./log/{data_name}/job_{suffix}.log
error = ./log/{data_name}/job_{suffix}.err
queue"""

def main():
    print data_name, txt_name
    out_dir = 'root://cms-xrdr.private.lo:2094///xrd/store/user/yyoun/DeepCMeson/analyser_'+data_name
    log_dir = './log/{}'.format(data_name)
    #Make out dir 
    if not os.path.isdir(log_dir):
        print 'Make log directory'
        os.mkdir(log_dir)
    else: 
        print "Folder is already exist: ",log_dir
        sys.exit()

    with open(txt_name, 'r') as txt_file:
        input_path_list = [each.rstrip('\n') for each in txt_file.readlines()]
    for in_path in input_path_list:
        in_name = os.path.basename(in_path)
        out_name = in_name.replace('delphes', 'out')
        suffix = in_name.lstrip("delphes_").rstrip(".root")
        out_path = os.path.join(out_dir, out_name)

        # write jds file
        with open('{}/submit.jds'.format(log_dir), 'w') as jds_file:
            jds_file.write(_JDS_FMT.format(in_path=in_path, out_path=out_path, suffix=suffix, data_name=data_name))
    
        command = "condor_submit -batch-name DeepCMeson_{} {}/submit.jds".format(data_name, log_dir)
        os.system(command)

if __name__ == '__main__':
    main()
