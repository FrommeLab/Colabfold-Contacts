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

                #combines file names and files together, removes linebreaks
                cat <(echo -n $outputn) $fileA <(echo -n ':') $fileB | tr -d '\n' | tr -d '*'> ./temp_file

                #adds header info needed for ColabFold
                #change output dir path (leave file name variable)!
                cat <(echo 'id,sequence') ./temp > /path/to/output/dir/$output
                
                #removes temp file, change output dir path
                rm ./temp_file

        done
                
done
