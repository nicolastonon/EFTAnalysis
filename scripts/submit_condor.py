#Adapted from Potato code
#[[USAGE EXAMPLE: ./scripts/submit_condor.py analyze]]

#!/usr/bin/env python
#from __future__ import print_function

#-----
# submit analyzer jobs
def analyze(args):
    print('== SUBMITTING ANALYZER JOB(S) ==\n')

    # process arguments
    test = args.test
    split = int(args.split)

    # create output directory
    from datetime import datetime
    outdir = 'd{}-t{}'.format(datetime.now().strftime('%Y%m%d'), datetime.now().strftime('%H%M%S'))
    output = 'condor/' + outdir

    import os
    if not os.path.exists('condor'):
        os.mkdir('condor')
    os.mkdir(output)

    #-- NB: may loop here to generate multiple jobs
    i = 0 #Default: submit 1 single job
    if  True:

        # construct job name
        jobname = 'test'

        # create jobinformation
        ji = {
            'name': jobname,
            'jobs': {},
        }
        i += 1

    # write jobinformation
    import json
    with open(output+'/jobs.json', 'w') as f:
        json.dump(ji, f, indent=4, separators=(',', ': '))

    # submit jobs
    if test:
        def call(text, **opts):
            print(text)
    else:
        from subprocess import call
    command = './scripts/job.sh {} ./analysis_main.exe'.format(os.getcwd())
    call('./scripts/cs.sh -d{} -n{} -a0:{} "{}"'.format(outdir, jobname, i-1, command), shell=True) #NB: use [-a0:i-1] to submit i jobs with arguments 0...i-1

    print('== DONE ! ==\n')
# end: submit analyzer jobs
#-----


#-----
# submit NN training jobs
def train(args):
    # process arguments
    test = args.test
    split = int(args.split)

    # create output directory
    from datetime import datetime
    outdir = 'd{}-t{}'.format(datetime.now().strftime('%Y%m%d'), datetime.now().strftime('%H%M%S'))
    output = 'condor/' + outdir
    import os
    if not os.path.exists('condor'):
        os.mkdir('condor')
    os.mkdir(output)

    #-- NB: may loop here to generate multiple jobs
    i = 0 #Default: submit 1 single job
    if  True:

        # construct job name
        jobname = 'test'

        # create jobinformation
        ji = {
            'name': jobname,
            'jobs': {},
        }
        i += 1

    # write jobinformation
    import json
    with open(output+'/jobs.json', 'w') as f:
        json.dump(ji, f, indent=4, separators=(',', ': '))

    # submit jobs
    from subprocess import call
    call('cp ./Train_Neural_Network.py {}'.format(output), shell=True) #Copy current executable to timestamped dir (so that it won't get modified anymore)

    pretend=''
    if test: pretend = ' --pretend' #cs.sh won't submit
    #command = '../scripts/job.sh {} python ./Train_Neural_Network.py'.format(os.getcwd())
    command = '../scripts/job.sh {} python ./'.format(os.getcwd())+output+'/Train_Neural_Network.py' #Execute timestamped, job-dependent executable
    print('EXECUTING: ', '../scripts/cs.sh -n{}{} -a0:{} "{}"'.format(jobname, pretend, i-1, command), '\n')
    call('../scripts/cs.sh -d{} -n{}{} -a0:{} "{}"'.format(outdir, jobname, pretend, i-1, command), shell=True) #NB: hardcoding relative path from NeuralNetworks. dir.
# end: submit training jobs
#-----

#-----
# check analyzer jobs
def check(args):
    # process arguments
    directories = args.directories
    no_resubmit = args.noResubmit
    keep_files = args.keepFiles

    # preparation of options: delete files
    import os
    if not keep_files:
        remove = os.remove
    else:
        def remove(name):
            print('rm', name)
    # preparation of options: submit jobs
    if not no_resubmit:
        from subprocess import call
    else:
        def call(name, **opts):
            print(name)

    # loop over all directories
    directories = filter(os.path.isdir, map(lambda d: d[:-1] if d[-1]=='/' else d, directories))
    for directory in directories:
        print('<<< Now at {}'.format(directory))

        # read job information
        infofile = '{}/jobs.json'.format(directory)
        if not os.path.exists(infofile):
            print('<<< Did not find jobinfo.json in {}, skip this directory!'.format(directory))
            continue
        import json
        with open(infofile) as f:
            jobinfo = json.load(f)
        modulename = jobinfo['module']
        jobname = jobinfo['name']
        jobs = jobinfo['jobs']
        analyze = jobinfo.get('analyze', None)

        # list of potential results files
        files = filter(lambda f: f.endswith('.root'), os.listdir(directory))

        # check for existing analyze results
        import ROOT
        from common.python.SampleMode import Mode, Sample
        if analyze:
            unfinished_analyze = []
            analyzername, filelistname, output, version = analyze
            for i in jobs:
                # check if output file of analyzer exists
                filename, sample, mode, start, skip = jobs[i]
                prefix = '{}-{}-{}-{}-job{}-'.format(
                    analyzername, version,
                    Sample.getName(sample), Mode.getName(mode),
                    i,
                )
                candidates = filter(lambda f: f.startswith(prefix), files)
                if not candidates:
                    unfinished_analyze.append(i)
                    continue
                if len(candidates)>1:
                    candidates.sort(reverse=True)
                    for c in candidates[1:]:
                        remove(directory+'/'+c)
                resultfile = directory+'/'+candidates[0]
                tfile = ROOT.TFile.Open(resultfile)
                passed = tfile and not tfile.IsZombie() and tfile.GetListOfKeys() and len(tfile.GetListOfKeys()) and not tfile.TestBit(ROOT.TFile.kRecovered)
                extrafiles = []
                if passed and tfile.Get('customResults'):
                    for key in tfile.Get('customResults').GetListOfKeys():
                        extrafiles.append(key.ReadObj().GetTitle())
                if tfile:
                    tfile.Close()
                extrafiles = [directory+'/'+ extrafile for extrafile in filter(lambda f: f.startswith("{0}-{1}".format(analyzername, version)) and "-job{0}-".format(i) in f, files)]
                extrafiles = filter(lambda f: f != resultfile, extrafiles)
                # print("result ",resultfile)
                # print("extra ",extrafiles)

                for extrafile in extrafiles:
                    tfile = ROOT.TFile.Open(extrafile)
                    passed = passed and tfile and not tfile.IsZombie() and tfile.GetListOfKeys() and len(tfile.GetListOfKeys()) and not tfile.TestBit(ROOT.TFile.kRecovered)
                    tfile.Close()
                if not passed:
                    unfinished_analyze.append(i)
                    remove(resultfile)
                    for extrafile in extrafiles:
                        remove(extrafile)
                    continue

            # print results of failed jobs
            if unfinished_analyze:
                print('<<< Analyzer jobs have died:', ', '.join(map(lambda i: '{}'.format(i), sorted(unfinished_analyze, key=int))))
            else:
                print('<<< All analyzer jobs are finished')

            # re-submit failed jobs
            command = './common/scripts/job.sh {} ./analysis_main.exe'.format(os.getcwd())
            for nJob in unfinished_analyze:
                call('common/scripts/cs.sh -n{0}-{2} "{1} {2}"'.format(
                    jobname, command, nJob
                ), shell=True)
# end: check analyzer jobs
#-----

#-----
# submit jobs via GUI
def GUI():
    print('opening GUI...')
    import common.python.SubmitterGui
# end: submit jobs via GUI
#-----

#-----
# parse command line arguments
def Submit():
    # setup argument parser
    import argparse
    import os, json
    parser = argparse.ArgumentParser(prog='./submit')
    subparsers = parser.add_subparsers(dest='subparser')

    # command line arguments: analyze
    parserA = subparsers.add_parser(
        'analyze',
        help='run an Analyzer on ntuples'
    )
    parserA.add_argument(
        '-t', '--test', action='store_true',
        help='don\'t submit, just test'
    )
    parserA.add_argument(
        '-s', '--split', default=1,
        help='split each job into multiple jobs'
    )

    # command line arguments: train
    parserT = subparsers.add_parser(
        'train',
        help='Train NN on ntuples'
    )
    parserT.add_argument(
        '-t', '--test', action='store_true',
        help='don\'t submit, just test'
    )
    parserT.add_argument(
        '-s', '--split', default=1,
        help='split each job into multiple jobs'
    )

    # command line arguments: check
    parserC = subparsers.add_parser(
        'check',
        help='check and resubmit Analyzer jobs to HTCondor'
    )
    parserC.add_argument(
        'directories', nargs='+', metavar='directory',
        help='directory that contains the jobinfo.json'
    )
    parserC.add_argument(
        '-n', '--noResubmit', action='store_true',
        help='don\'t resubmit jobs'
    )
    parserC.add_argument(
        '-k', '--keepFiles', action='store_true',
        help='don\'t delete any files'
    )

    # command line arguments: GUI
    parserGUI = subparsers.add_parser(
        'GUI',
        help='open primitive GUI'
    )

    args = parser.parse_args()


    # parse arguments and call subparser
    if args.subparser == 'analyze':
        analyze(args)
    elif args.subparser == 'train':
        train(args)
    elif args.subparser == 'check':
        check(args)
    elif args.subparser == 'GUI':
        GUI()

if __name__ == '__main__':
    Submit()
# end: parse command line arguments
#-----
