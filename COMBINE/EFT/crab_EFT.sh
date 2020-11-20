set -x
set -e
ulimit -s unlimited
ulimit -c 0

function error_exit
{
  if [ $1 -ne 0 ]; then
    echo "Error with exit code ${1}"
    if [ -e FrameworkJobReport.xml ]
    then
      cat << EOF > FrameworkJobReport.xml.tmp
      <FrameworkJobReport>
      <FrameworkError ExitStatus="${1}" Type="" >
      Error with exit code ${1}
      </FrameworkError>
EOF
      tail -n+2 FrameworkJobReport.xml >> FrameworkJobReport.xml.tmp
      mv FrameworkJobReport.xml.tmp FrameworkJobReport.xml
    else
      cat << EOF > FrameworkJobReport.xml
      <FrameworkJobReport>
      <FrameworkError ExitStatus="${1}" Type="" >
      Error with exit code ${1}
      </FrameworkError>
      </FrameworkJobReport>
EOF
    fi
    exit 0
  fi
}

trap 'error_exit $?' ERR

if [ $1 -eq 1 ]; then
  ./combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt -t -1 -v 2  -M MultiDimFit -d EFTWorkspace.root --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-8,8:cpt=-8,8 --points 1000 --firstPoint 0 --lastPoint 999 -n .EFT.POINTS.0.999
fi
tar -cf combine_output.tar higgsCombine*.root
rm higgsCombine*.root

