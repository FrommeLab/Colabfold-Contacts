# Colabfold-Pipeline

`parse_fastas.py` takes an input fasta file containing bait and target sequences and writes out files to the local directory that are ready for 
local colabfold complex prediction. Requires biopython and expects one of the sequences to be named >bait
<pre>python parse_fastas.py --i inputs.fasta</pre>
