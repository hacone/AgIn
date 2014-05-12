##Installation
    $ mkdir AgIn
    $ cd AgIn
    $ git clone git://github.com/hacone/AgIn.git
    $ ./install.sh
You can run the demo by `$ ./run.sh HNI_Tol2/predict`

##Analyse your data
To run the analysis for your data, you must prepare folloing input
files and specify your task file.

####1. Place reference files (`.fasta` and `.fasta.fai`) at `input/`
Then your reference file must be specified in `Reference` field of
`input/config.dat` in relative path from `input/`
*Hint: When `Reference` field is missing, the default behaviour is to
regard `input/GenRef.fasta(.fai)` as symbolic links to your
reference.*

####2. Place IPD files at `input/IPD/`
IPD files can be prepared by processing `modification.csv` file
produced by `P_ModificationDetection` module of SMRT Pipe.
Specifically, an IPD file is CSV containing, for each position,
folloing columns:
    `Strand(0 or 1), Template position, mean IPD ratio, SMRT read coverage`
The columns must be in this order; the rows are to be sorted in 1st,
and then 2nd column (cf. our samples in `input/IPD/`).
Each IPD file must represent a single reference entry,
e.g., you will have 26 IPD files for human methylome analysis.

####3. Create a task file as `task/<YourTaskName>.json`
To do this, please copy and modify our sample task file
`task/HNI_Tol2/predict.json`.
You have to specify input files and desired prefix of output files.
Input files are specified by 2 elements list of reference entry name
(as in `.fai` file) and corresponding IPD file (relative path from
`input/IPD/`).
Output prefix would be anything you like and may contain `/` if you
want to get outputs in a specific directory.

####4. Set the parameters (optional)
If you want to change the parameters for prediction, please edit
`input/config.dat`.
E.g., by setting:

    L = 100
    Gamma = -1.70

AgIn uses `100` as minimum number of CpG sites of predicted regions
and `-1.70` as decision threshold for prediction.

####5. Run
You can now run your analysis by:
`$ ./run.sh <YourTaskName>`

## Output files specification
Names of output files would look like:

    <prefix><reference entry name>.allpred.dat
    <prefix><reference entry name>.avgscr.dat

`.allpred.dat` contains 4 columns:

    1: position of CpG site
    2: its raw prediction score
    3: predicted class(0 for hypo, 1 for hyper)
    4: average SMRT coverage of 21 bp window

If you are interested only in 1st and 3rd columns of these,
each predicted region is reported in more conprehensive `.avgscr.dat`
which also contains 4 columns:

    1: position of first CpG site
    2: position of last CpG site
    3: average score
    4: number of CpG sites contained in the region

## Lisence

Copyright of this work `AgIn` is held by Yuta Suzuki
(ysuzuki@cb.k.u-tokyo.ac.jp) and  currently distributed under  ?????
license. (NOTE: subject to change)

Copyright of `src/script/launch` is held by Taro L.
Saito(leo@xerial.org) and licensed under Apache License Version 2.0
(http://www.apache.org/licenses/LICENSE-2.0)
