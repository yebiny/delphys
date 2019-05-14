import os

#JobName = DeepC

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
output = TestJob/job_{suffix}.log
error = TestJob/job_{suffix}.err
queue"""


def main():
    out_dir = 'root://cms-xrdr.private.lo:2094///xrd/store/user/yyoun/xrootd/DeepCMeson/2-Analysis.v2'

    with open("dataset_DeepC.txt", 'r') as txt_file:
        input_path_list = [each.rstrip('\n') for each in txt_file.readlines()]

    for in_path in input_path_list:
        in_name = os.path.basename(in_path)
        out_name = in_name.replace('delphes', 'out')

        suffix = in_name.lstrip("delphes_").rstrip(".root")

        out_path = os.path.join(out_dir, out_name)

        # write jds file
        with open('DatasetV2/submit.jds', 'w') as jds_file:
            jds_file.write(_JDS_FMT.format(in_path=in_path, out_path=out_path, suffix=suffix))
            
        command = "condor_submit -batch-name DeepCMeson DatasetV2/submit.jds"


        os.system(command)

if __name__ == '__main__':
    main()
