import rdkit.Chem.Scaffolds.MurckoScaffold
import scaffoldgraph as sg
from rdkit import Chem


def mol_equal(s, t):
    return s.HasSubstructMatch(t) and t.HasSubstructMatch(s)


def compute_scaffold(s):
    s = Chem.MolFromSmiles(s)
    scaff = Chem.Scaffolds.MurckoScaffold.GetScaffoldForMol(s)
    return [Chem.MolToSmiles(scaff)]


def compute_lower(s):
    sorig = s
    s = Chem.MolFromSmiles(s)
    results = map(lambda x: Chem.MolToSmiles(x), sg.get_next_murcko_fragments(s))
    results = list(results)
    check_lower(sorig, results[0])
    return list(results)


def check_scaffold(s, t):
    s = Chem.MolFromSmiles(s)
    t = Chem.MolFromSmiles(t)

    if s is None or t is None:
        return -1

    try:
        scaffold = Chem.Scaffolds.MurckoScaffold.GetScaffoldForMol(s)
    except:
        return -1

    if mol_equal(scaffold, t):
        return 1
    else:
        return 0


def check_lower(s, t):
    s = Chem.MolFromSmiles(s)
    t = Chem.MolFromSmiles(t)
    if s is None or t is None:
        return -1
    if s.HasSubstructMatch(t):
        return 1
    else:
        return 0


def check_upper(s, t):
    s = Chem.MolFromSmiles(s)
    t = Chem.MolFromSmiles(t)

    if s is None or t is None:
        return -1
    if t.HasSubstructMatch(s):
        return 1
    else:
        return 0


def check_expand(s, t):
    s = Chem.MolFromSmiles(s)
    t = Chem.MolFromSmiles(t)

    if s is None or t is None:
        return -1
    try:
        scaffold = Chem.Scaffolds.MurckoScaffold.GetScaffoldForMol(t)
    except:
        return -1

    if mol_equal(scaffold, s):
        return 1
    else:
        return 0
