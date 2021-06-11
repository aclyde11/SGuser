## ScaffoldLattice REPL

This is sample code to provide an interactive environment which interfaces a transformer model (as a compression of
chemical space).

## Setup
The provided env file should get you started. If problems encountered, just install packages as needed
with python.
```shell
conda env create -f environment.yml #create conda env for SGUser
conda activate SGUser
python sg_repl.py out.xlsx 
```

## Usage
Scaffold graph as a very simple syntax which uses reverse polish notation to evaluation expressions. 
Expressions have a simple syntax
```
CMDS = [LOWER, UPPER, EXPAND, SCAFFOLD]
CMD SMILES -> SMILES
SMILES -> SMILES
```
An exciting part of this project is that depending on the command used, the code running in the background
may be a CPU algorithm, or it may be a call to a generative transformer model. 

Here is an example command to get some molecules which contain the scaffold of the following compound
```shell
python SGuser/sg_repl.py out.xlsx
ScaffoldLattice REPL
To run with history file (which includes pictures) please run python sg_repl.py [save file .html or .xlsx]
	The syntax is 'CMD2 CMD1 SMILES' where CMDN is LOWER,UPPER,EXPAND,SCAFFOLD
	To exit, type exit
	To change the sampling type 'SETN N' where N is the number of samples to draw
 from an underlying set produced by UPPER, LOWER, and EXPAND
Version 0.1 
contact aclyde@uchicago.edu
>>> EXPAND UPPER SCAFFOLD O=C(Oc1cncc(Cl)c1)c1cccc2[nH]ccc12
	O=C(Oc1cccnc1)c1cccc2[nH]ccc12
INFO: EXECUTED SCAFFOLD

INFO: NOW EXECUTING UPPER on RESULTS[0] O=C(Oc1cccnc1)c1cccc2[nH]ccc12

		O=C(Oc1cccnc1)c1cccc2c1ccn2Cc1ccccc1
		O=C(Oc1cccnc1)c1cccc2c1ccn2CC1CC1
		O=C(Oc1cccnc1)c1cccc2c1ccn2CC#Cc1ccccc1
		O=C(Oc1cccnc1)c1cccc2c1ccn2CC1CCNCC1
		O=C(Oc1cccnc1)c1cccc2[nH]c3ccccc3c12
INFO: EXECUTED UPPER

INFO: NOW EXECUTING EXPAND on RESULTS[0] O=C(Oc1cccnc1)c1cccc2c1ccn2Cc1ccccc1

			Cc1cccc(Cn2cc(C)c3c(C(=O)Oc4cccnc4)cccc32)c1
			Cc1ccc(Cn2ccc3c(C(=O)Oc4cccnc4)cccc32)cc1
			Cc1cc(C)cc(Cn2ccc3c(C(=O)Oc4cccnc4)c(Cl)ccc32)c1
			COc1ncccc1OC(=O)c1cccc2c1ccn2Cc1ccc(C)cc1
			Cc1cc(C(=O)Oc2cccnc2)c2ccn(Cc3cccc(Cl)c3Br)c2c1
INFO: EXECUTED EXPAND
```
You will notice that each command is executed separately. This is because most commands
besides SCAFFOLD generate a set lazily. In order to execute another command
that set must be sampled to a particular compound to run on the next query.
By default, and currently is not changeable, the first result of the prior run is
used to continue the chain of execution. At the enb of running, you will get a nice
output document which contains pictures of everything you did.



