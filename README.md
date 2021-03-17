## PredGPI: prediction of GPI-anchors in proteins based on HMMs and SVMs

#### Publication
Pierleoni A, Martelli PL, Casadio R. PredGPI: a GPI-anchor predictor. BMC
Bioinformatics. 2008 Sep 23;9:392.

#### Installation and configuration

Before running predgpi you need to set and export a variable named PREDGPI_HOME to point to the program installation dir:
```
$ export PREDGPI_HOME='/path/to/predgpi'
```

#### Usage

To run the program:

```
./predgpi.py -f FASTA -o OUTF -m {json|gff3}

```
