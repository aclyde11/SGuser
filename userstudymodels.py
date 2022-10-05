from rdkit import Chem
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem import rdqueries
from rdkit.Chem.rdMolDescriptors import CalcNumHeteroatoms
from numpy import random
from rdkit.Chem.Lipinski import NumHDonors

def model1(smiles):
    mol = Chem.MolFromSmiles(smiles)
    q = rdqueries.AtomNumEqualsQueryAtom(6)
    NC = len(mol.GetAtomsMatchingQuery(q))
    NHET = CalcNumHeteroatoms(mol)
    return 1.46+ 0.11 * NC - 0.11 * NHET

def model2(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return NumHDonors(mol) * random.rand()

def model3(smiles):
    return logp(smiles) * 2*random.rand()

def logp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return MolLogP(mol)

if __name__ == '__main__':
    txt = "COc1cc(NC(=O)CN[C@@H]2CC[C@@H](C)[C@H](C)C2)cc(OC)c1"
    print(model1(txt))
    print(model2(txt))
    print(model3(txt))
    print(logp(txt))