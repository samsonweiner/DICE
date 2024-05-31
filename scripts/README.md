# Evaluation scripts

The `execute.py` script calls functions located in `utilities.py` to execute the methods used in the DICE evaluations. Running the scripts to evaluate external methods requires they are installed. It also requires the following python packages: ete3, skbio. The scripts take as input the copy number profiles in the same format as DICE, and automatically applies transformations so the data is in the correct format to run each method.

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
