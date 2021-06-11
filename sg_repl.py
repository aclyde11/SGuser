import os
import sys

os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'  # some weird libomp.dylib error fix

import subprocess
import tempfile
import selfies
from rdkit import Chem
from rdkit.Chem import AllChem, PandasTools
import re

from rdkit import RDLogger
import pandas as pd

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
RDLogger.DisableLog('rdApp.info')

MODEL = "saved_models/algebraDataSmall_step_199000.pt"
N = 6
HISTORY = []

notnone = lambda x: x is not None

import validate as molopsval


def check_cmd(op, inputsmi, results):
    checked = []
    for result in results:
        if op == 'SCAFFOLD':
            if molopsval.check_scaffold(inputsmi, result) == 1:
                checked.append(result)
        elif op == 'EXPAND':
            if molopsval.check_expand(inputsmi, result) == 1:
                checked.append(result)
        elif op == 'LOWER':
            if molopsval.check_lower(inputsmi, result) == 1:
                checked.append(result)
        elif op == 'UPPER':
            if molopsval.check_upper(inputsmi, result) == 1:
                checked.append(result)
    return checked


def tokenzie_smile(smi):
    pattern = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    assert smi == ''.join(tokens)
    return ' '.join(tokens)


def tokenzie_selfie(smi):
    tokens = selfies.split_selfies(smi)
    return ' '.join(tokens)


def decode_selfies(selfie):
    results = []
    for s in selfie:
        try:
            results.append(selfies.decoder(s))
        except Exception as e:
            print("Error ", e, s)
    return results


def make_call(cmd, useselfie=False):
    with tempfile.TemporaryDirectory() as tempdir:
        with open(f"{tempdir}/t1.txt", 'w') as f:
            op, smiles = cmd.strip().split(' ')
            if useselfie:
                selfie = selfies.encoder(smiles)
                selfie = tokenzie_selfie(selfie)
                outcmd = f"{op} {selfie}\n"
            else:
                smiles_toks = tokenzie_smile(smiles)
                outcmd = f"{op} {smiles_toks}\n"
            for i in range(N):
                f.write(outcmd)
        if op in ['UPPER', 'EXPAND']:
            cmdstr = f"onmt_translate -model {MODEL} -src {tempdir}/t1.txt -output {tempdir}/t2.txt --batch_size 4" \
                     f" --random_sampling_topk -1 --beam_size 5".split(" ")
            subprocess.run(cmdstr, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

            with open(f"{tempdir}/t2.txt", 'r') as f:
                results = list(map(lambda x: x.strip().replace(' ', ''), f))
                if useselfie:
                    results = decode_selfies(results)
                results_valid = set()
                for result in results:
                    m = Chem.MolFromSmiles(result)
                    if m is not None:
                        results_valid.add(Chem.MolToSmiles(m))
        else:
            if op == 'LOWER':
                results_valid = molopsval.compute_lower(smiles)
            elif op == 'SCAFFOLD':
                results_valid = molopsval.compute_scaffold(smiles)
    return list(results_valid)[:N]


def info():
    print("ScaffoldLattice REPL")
    print("To run with history file (which includes pictures) please run python sg_repl.py [save file .html or .xlsx]")
    print("\tThe syntax is '>>> CMD2 CMD1 SMILES' where CMDN is LOWER,UPPER,EXPAND,SCAFFOLD")
    print("\tTo exit, type exit")
    print(
        "\tTo change the sampling type 'SETN N' where N is the number of samples to draw\n from an underlying set produced by UPPER, LOWER, and EXPAND")
    print("Version 0.1 \ncontact aclyde@uchicago.edu")


def checkvalid(smiles):
    try:
        n = Chem.MolFromSmiles(smiles)
        return (n is not None)
    except:
        return False


def gethash(smiles):
    hashes = set()
    keepers = []
    for smile in smiles:
        try:
            n = Chem.MolFromSmiles(smile)
            if n is None:
                continue
            fp1 = AllChem.GetMorganFingerprint(n, 2)
            if fp1 in hashes:
                continue
            keepers.append(smile)
            hashes.add(fp1)
        except:
            continue
    return keepers


def parsecmd(x):
    original_command = x
    x = x.strip()
    if "SETN" in x:
        N = int(x.split(" ")[-1])
        return [f"Set sampling n to {N}."]
    smiles = x.split(" ")[-1]
    x = x.split(" ")[:-1]

    for idx, op in enumerate(x[::-1]):
        resulting_command = f"{op} {smiles}"
        results = make_call(resulting_command)
        # print("Results from call", results)
        results = gethash(results)
        # print("Unqiue and valid", results)
        results = check_cmd(op, smiles, results)
        # print("correct", results)
        for result in results:
            history_node = {'Original Command': x, 'Executed Command': resulting_command, 'SMILES': result}
            print("".join((idx + 1) * ['\t']) + f"{result}")
            HISTORY.append(history_node)
        if len(x) > 1:
            smiles = results[0]
            print(f"INFO: EXECUTED {x[len(x) - 1 - idx]}\n")
            if len(x) - 2 - idx >= 0:
                print(f"INFO: NOW EXECUTING {x[len(x) - 2 - idx]} on RESULTS[0] {smiles}\n")


def repl(save_file=None):
    info()
    try:
        while True:
            try:
                _in = input(">>> ")
                if "exit" in _in:
                    if save_file is not None:
                        history = pd.DataFrame.from_dict(HISTORY)
                        print(history)
                        PandasTools.AddMoleculeColumnToFrame(history, 'SMILES', molCol='Molecule')
                        PandasTools.ChangeMoleculeRendering(history, 'PNG')
                        if save_file.split(".")[-1] == 'html':
                            htm = history.to_html()
                            with open(save_file, 'w') as f:
                                f.write(htm)
                        elif save_file.split(".")[-1] == 'xlsx':
                            PandasTools.SaveXlsxFromFrame(history, save_file, molCol='Molecule')
                        else:
                            print(f"Provided save file {save_file} however this file must end in .html or .xlsx.")
                    exit()
                try:
                    parsecmd(_in)
                except KeyboardInterrupt as e:
                    print("\nExiting...")
            except KeyboardInterrupt:
                print("Press Ctrl-C to terminate while statement")
                break
    except KeyboardInterrupt as e:
        print("\nExiting...")
    print("Exiting SG REPL.")


if __name__ == '__main__':
    if len(sys.argv) > 1:
        repl(sys.argv[1])
    else:
        repl()
