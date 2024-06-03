# Evaluation scripts

Scripts are provided to execute the methods used in the DICE evaluation, including [DICE](https://github.com/samsonweiner/DICE), [MEDICC2](https://bitbucket.org/schwarzlab/medicc2/src), [MEDALT](https://github.com/KChen-lab/MEDALT), [cnp2cnp](https://github.com/AEVO-lab/cnp2cnp), [sitka](https://github.com/UBC-Stat-ML/sitkatree), and [lazac](https://github.com/raphael-group/lazac-copy-number). The `execute.py` script calls functions located in `utilities.py` to run these methods. Using the scripts requires the following additional python packages: ete3, skbio. Additionally, the external methods must be installed. The path to each installed method can be specified at the top of the `execute.py` script. It takes as input the copy number profiles in the same format as DICE, and automatically applies transformations to the input data to match the required format necessary for running each method.

The options are as follows:

```
    -i, --input: Path to input copy number profiles.
    -p, --prepare-input: Convert the input data to the correct formats to run each method.
    -n, --num-processor: Number of cores to use in parallel (if methods allow). 
    -d, --dice: Run all DICE variants.
    -0, --medicc2: Run MEDICC2.
    -1, --medalt: Run MEDALT.
    -2, --cnp2cnp: Run cnp2cnp.
    -3, --sitka: Run sitka.
    -4, --lazac: Run lazac.
```

