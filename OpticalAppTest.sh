#!/bin/bash -l 
usage(){ cat << EOU
~/OpticalApp/OpticalAppTest.sh
================================

Performs an optical simulation and dumps optical photon
velocities within photon history categories. 

::

    VERS=ENV  ~/OpticalApp/OpticalAppTest.sh
    VERS=1042 ~/OpticalApp/OpticalAppTest.sh
    VERS=1120 ~/OpticalApp/OpticalAppTest.sh

    BP=G4OpBoundaryProcess::PostStepDoIt VERS=1042 ~/OpticalApp/OpticalAppTest.sh
    ## find dbg only works with VERS=ENV or VERS=1042 where 1042 matches the default version 

    BP=G4OpAbsorption::PostStepDoIt ~/OpticalApp/OpticalAppTest.sh
    BP=G4OpAbsorption::GetMeanFreePath ~/OpticalApp/OpticalAppTest.sh


Started from Opticks example ~/o/examples/Geant4/OpticalApp/OpticalAppTest.sh 


* TODO: more realistic physics : absorption, reemission


EOU
}

cd $(dirname $(realpath $BASH_SOURCE))

name=OpticalAppTest
export FOLD=/tmp/$name
mkdir -p $FOLD

bin=$FOLD/$name
script=$name.py 


gdb__ () 
{ 
    : opticks/opticks.bash prepares and invokes gdb - sets up breakpoints based on BP envvar containing space delimited symbols;
    if [ -z "$BP" ]; then
        H="";
        B="";
        T="-ex r";
    else
        H="-ex \"set breakpoint pending on\"";
        B="";
        for bp in $BP;
        do
            B="$B -ex \"break $bp\" ";
        done;
        T="-ex \"info break\" -ex r";
    fi;
    local runline="gdb $H $B $T --args $* ";
    echo $runline;
    date;
    eval $runline;
    date
}


unset JUNOTOP

vers=ENV
VERS=${VERS:-$vers}

if [ "$VERS" == "ENV" ]; then 

    echo $BASH_SOURCE : use geant4_config from environment 

elif [ "$VERS" == "1042" -o "$VERS" == "1120" ]; then

    if [ "$VERS" == "1120" ]; then 

        base=/cvmfs/juno.ihep.ac.cn/centos7_amd64_gcc1120/Pre-Release/J24.1.x-g411/ExternalLibs
        my_clhep_prefix=$base/CLHEP/2.4.7.1
        my_geant4_prefix=$base/Geant4/11.2.0

    elif [ "$VERS" == "1042" ]; then 

        base=/cvmfs/juno.ihep.ac.cn/centos7_amd64_gcc1120/Pre-Release/J23.1.0-rc3/ExternalLibs
        my_clhep_prefix=$base/CLHEP/2.4.1.0
        my_geant4_prefix=$base/Geant4/10.04.p02.juno
    fi

    source $my_clhep_prefix/bashrc
    source $my_geant4_prefix/bashrc
fi 

G4CFG=$(which geant4-config)

vars="BASH_SOURCE name bin FOLD G4CFG"


#export OpticalApp__GeneratePrimaries_DEBUG_GENIDX=50000
#export OpticalApp__PreUserTrackingAction_UseGivenVelocity_KLUDGE=1 

#export OpticalApp__InitRunManager_OpticalPhysics=G4OpAbsorption,G4OpBoundaryProcess
export OpticalApp__InitRunManager_OpticalPhysics=G4OpBoundaryProcess,G4OpAbsorption



# -Wno-deprecated-copy \

defarg=info_build_run
[ -n "$BP" ] && defarg=info_build_dbg 

arg=${1:-$defarg}

if [ "${arg/info}" != "$arg" ]; then
   for var in $vars ; do printf "%20s : %s \n" "$var" "${!var}" ; done 
fi


opt=""
case $(uname) in 
  Linux) opt="-Wno-deprecated-copy" ;;
  Darwin) opt="-Wno-deprecated" ;;
esac


if [ "${arg/build}" != "$arg" ]; then
    gcc $name.cc \
            -I. \
            -g \
            $(geant4-config --cflags) \
            -Wno-shadow $opt \
            $(geant4-config --libs) \
            -lstdc++ -lm \
            -o $bin
    [ $? -ne 0 ] && echo $BASH_SOURCE : build error && exit 1 
fi

if [ "${arg/run}" != "$arg" ]; then
    $bin
    [ $? -ne 0 ] && echo $BASH_SOURCE : run error && exit 2
fi

if [ "${arg/dbg}" != "$arg" ]; then
   
    case $(uname) in 
       Linux) gdb__ $bin ;;
       Darwin) lldb__ $bin ;;
    esac
    [ $? -ne 0 ] && echo $BASH_SOURCE : dbg error && exit 3
fi

exit 0 

