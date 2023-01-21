#!/bin/sh

phits_MPI="/home/dnw/phits327A/phits/bin2/phits_LinGfort_MPI"
phits_OMP="/home/dnw/phits327A/phits/bin2/phits_LinGfort_OMP"
phits_single="/home/dnw/phits327A/phits/bin2/phits_LinGfort"
mpi_command="mpirun"


while getopts "M:Ohs:" option; do
    case "${option}" in
	\?) echo "See -h for help."
	    exit;;

	M)  phitsin_path="$(dirname -- "$(readlink -f -- "$phits_MPI")")/phits.in"
	    nextopt=${OPTARG}
	    if [[ $nextopt != ${@: -1} ]] ; then
		np=$nextopt
	    else
		np=$(nproc --all)
	    fi
	    echo "file = $(readlink -f ${@: -1})" > ${phitsin_path}
	    eval "${mpi_command} -np ${np} ${phits_MPI}"

	    exit;;

	O) ${phits_OMP} < ${@: -1}
	   exit;;

	h) echo "phits: nicer command-line interface to JAEA's Particle Heavy Ion Transport code System"
	   echo "Usage: phits -[MOhs] file"
	   echo "Options:"
	   echo "    none: run the standard PHITS single-core binary on file"
	   echo "    -M [proc]: copy the path of file to phits.in and run a PHITS MPI binary on all cores of the current machine. If specified, use [proc] processes; otherwise, use every core."
	   echo "    -O: run a PHITS OpenMP binary on file."
	   echo "    -h: display this message."
	   echo "    -s ['omp'|'mpi'|'single']: switch the PHITS binary of the given type run by this command to that specified in file."

	   exit;;

	s) s=${OPTARG}
	   file=${@: -1}
	   if [[ $s = $file ]] ; then
	       echo "No new binary specified for option -s."
	   elif [[ $s = 'mpi' ]] ; then
	       sed -i "s/^phits_MPI=/phits_MPI=${file}" $BASH_SOURCE
	   elif [[ $s = 'omp' ]] ; then
	       sed -i "s/^phits_OMP=/phits_OMP=${file}" $BASH_SOURCE
	   elif [[ $s = 'single' ]] ; then
	       sed -i "s/^phits_single=/phits_single=${file}" $BASH_SOURCE
	   else
	       echo "Unrecognized binary type \"${s}\"."
	   fi

	   exit;;
    esac
done

shift $(($OPTIND - 1))

# no option given; use single-thread binary
if [[ -n $1 ]] ; then
    ${phits_single} < $1
else
    echo "Error: no input file specified."
fi
