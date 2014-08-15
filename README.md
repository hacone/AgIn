    *NOTE: These instructions are subject to change*

##Installation
    $ git clone git://github.com/hacone/AgIn.git AgIn

##Analyse your data
You need a `modification.csv` file, which is produced by `P_ModificationDetection` module of SMRT Pipe.
Then run as:

    $ target/dist/bin/launch -i path/to/modifications.csv -f path/to/reference.fasta -o output_prefix -g gamma -l min_length predict

Then you will obtain three output files:

    output_prefix.gff
    output_prefix_class.wig
    output_prefix_coverage.wig

Please make sure the index file for the specified fasta file be placed in the same directory as fasta.


####The parameters
If you want to change the parameters for prediction, supply them in the command-line options;
By you specifying `-l 100 -g -1.70`, 
then  AgIn uses `100` as minimum number of CpG sites the predicted regions must contain,
and `-1.70` as decision threshold for prediction.

## Output files specification
Though these are fairly intuitive, some description would be added here soon.

## Lisence

Copyright of this work `AgIn` is held by Yuta Suzuki
(ysuzuki@cb.k.u-tokyo.ac.jp) and currently distributed under the (modified) BSD license.

Copyright of `src/script/launch` and some other scripts is held by Taro L.
Saito(leo@xerial.org) and licensed under Apache License Version 2.0
(http://www.apache.org/licenses/LICENSE-2.0)
