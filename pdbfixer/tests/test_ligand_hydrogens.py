"""Tests for PDBFixer.addHydrogensToLigands (RDKit, template-driven, offline)."""
import io

import pytest

pytest.importorskip("rdkit")

from rdkit import Chem
from rdkit.Chem import AllChem
from openmm import app, unit, Vec3

import pdbfixer


# (name, SMILES): a spread of chemistries -- ester+acid, fused aromatics with
# carbonyls, a large drug, an amidine, and a halogenated amide+sulfone.
LIGANDS = [
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("imatinib", "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"),
    ("benzamidine", "c1ccc(cc1)C(=N)N"),
    ("halo_amide", "O=C(N)c1cc(F)ccc1OCC1CCN(C(=O)C)CC1"),
]


def _ref_h_count(smi):
    base = Chem.MolFromSmiles(smi)
    return Chem.AddHs(base).GetNumAtoms() - base.GetNumAtoms()


def _fixer_from_heavy_atoms(smi, resname="LIG", seed=7):
    """Embed a molecule, strip its hydrogens, and expose only the heavy atoms
    (coordinates + elements) through PDBFixer -- mimicking a generated ligand
    written to a PDB file with no hydrogens."""
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(m, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(m)
    heavy = Chem.RemoveHs(m)
    conf = heavy.GetConformer()

    top = app.Topology()
    res = top.addResidue(resname, top.addChain())
    coords = []
    for i, atom in enumerate(heavy.GetAtoms()):
        el = app.element.Element.getByAtomicNumber(atom.GetAtomicNum())
        top.addAtom(f"{atom.GetSymbol()}{i}", el, res)  # unique names
        p = conf.GetAtomPosition(atom.GetIdx())
        coords.append(Vec3(p.x * 0.1, p.y * 0.1, p.z * 0.1))  # Angstrom -> nm

    buf = io.StringIO()
    app.PDBFile.writeFile(top, coords * unit.nanometer, buf)
    buf.seek(0)
    return pdbfixer.PDBFixer(pdbfile=buf)


@pytest.mark.parametrize("name,smi", LIGANDS)
def test_added_count_matches_reference(name, smi):
    """Each ligand gains exactly the hydrogens its SMILES implies, with a
    consistent topology and the hydrogens actually present."""
    fixer = _fixer_from_heavy_atoms(smi)
    added = fixer.addHydrogensToLigands({"LIG": smi})
    ref = _ref_h_count(smi)
    topo_h = sum(1 for a in fixer.topology.atoms()
                 if a.element is not None and a.element.symbol == "H")
    assert added.get("LIG", 0) == ref
    assert topo_h == ref
    assert fixer.topology.getNumAtoms() == len(fixer.positions)


def test_accepts_rdkit_mol_template():
    """The template may be an RDKit Mol, not only a SMILES string."""
    name, smi = LIGANDS[0]
    fixer = _fixer_from_heavy_atoms(smi)
    added = fixer.addHydrogensToLigands({"LIG": Chem.MolFromSmiles(smi)})
    assert added.get("LIG", 0) == _ref_h_count(smi)


def test_no_template_leaves_ligand_untouched():
    """A ligand without a matching template keeps zero hydrogens."""
    smi = LIGANDS[1][1]
    fixer = _fixer_from_heavy_atoms(smi)
    before = fixer.topology.getNumAtoms()
    assert fixer.addHydrogensToLigands({}) == {}
    assert fixer.addHydrogensToLigands({"OTHER": smi}) == {}
    assert fixer.topology.getNumAtoms() == before


def test_mismatched_template_is_rejected():
    """A template with the wrong heavy-atom count is reported, not applied."""
    smi = LIGANDS[0][1]
    fixer = _fixer_from_heavy_atoms(smi)
    before = fixer.topology.getNumAtoms()
    added = fixer.addHydrogensToLigands({"LIG": "CCO"})  # 3 heavy atoms, wrong
    assert added.get("LIG", 0) == 0
    assert fixer.topology.getNumAtoms() == before
