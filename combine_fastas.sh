#!/bin/bash

#RCV
#08.30.22

#combines fasta files in dirA and dirB into Colabfold format

#input directories
#change to desired folders (e.g., all SNAREs in dirA and all other proteins in dirB)
dirA=/path/to/inputsA/
dirB=/path/to/inputsB/

#loops through all files in dirA
for fileA in $dirA/*.fas
        do

#loops through all files in dirB
        for fileB in $dirB/*.fas
                do

                #assings variables for names, janky I know
                extA=$(basename -- "$fileA")
                outputA="${extA%.*}"
                extB=$(basename -- "$fileB")
                outputB="${extB%.*}"
                outputn="$outputA"_"$outputB"","
                output="$outputA"_"$outputB"".csv"

                #combines file names and files together, removes linebreaks, change output dir path (leave file name as temp)
                cat <(echo -n $outputn) $fileA <(echo -n ':') $fileB | tr -d '\n' | tr -d '*'> /path/to/output/dir/temp

                #adds header info needed for ColabFold, change output dir path (leave file name as temp)
                cat <(echo 'id,sequence') /path/to/output/dir/temp > /home/rcv42/all_ORFs/golgiVsSnares/$output
                
                #removes temp file, change output dir path
                rm /path/to/output/dir/temp

        done
                
done
