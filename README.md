[![PDBFixer Continuous Integration Workflow](https://github.com/openmm/pdbfixer/actions/workflows/CI.yml/badge.svg)](https://github.com/openmm/pdbfixer/actions/workflows/CI.yml)

PDBFixer
========

PDBFixer is an easy to use application for fixing problems in Protein Data Bank files in preparation for simulating them.  It can automatically fix the following problems:

- Add missing heavy atoms.
- Add missing hydrogen atoms (proteins, chain caps, and — in this fork — arbitrary ligands).
- Build missing loops.
- Convert non-standard residues to their standard equivalents.
- Select a single position for atoms with multiple alternate positions listed.
- Delete unwanted chains from the model.
- Delete unwanted heterogens (with independent control over water and ligands).
- Build a water box for explicit solvent simulations.
- Energy-minimize the prepared structure.

See our [manual](https://htmlpreview.github.io/?https://github.com/openmm/pdbfixer/blob/master/Manual.html)

## Fork additions

This fork (QuantaBricks/pdbfixer) adds a few capabilities on top of upstream PDBFixer.
All existing method signatures are unchanged; the items below are additive (the only
removed method relative to earlier fork revisions is `addMissingHydrogens2`).

### Hydrogenating ligands with RDKit — `addHydrogensToLigands`

`addMissingHydrogens` protonates standard polymer residues **and protein chain caps**
(`ACE`, `NME`, `NH2`, `NMA`, `FOR`), but deliberately leaves ligands/heterogens
untouched: the hydrogen topology of an arbitrary small molecule cannot be recovered
reliably from heavy-atom coordinates alone.

`addHydrogensToLigands` fills that gap for molecules whose chemistry you already know
(e.g. ligands you generated yourself). You supply a reference structure per residue
name; RDKit transfers its bond orders onto the matching heavy atoms and places
hydrogens at 3D positions. It runs **fully offline** — no Chemical Component Dictionary
or network access is used.

```python
from pdbfixer import PDBFixer

fixer = PDBFixer(filename='complex.pdb')
fixer.addMissingHydrogens(7.0)                          # protein + chain caps
fixer.addHydrogensToLigands({'LIG': 'CC(=O)Oc1ccccc1C(=O)O'})  # ligand, via RDKit
```

The template value may be a SMILES string or an RDKit `Mol`, and must contain exactly
the same heavy atoms as the residue in the model (hydrogens in the reference are
ignored). A ligand with no matching template is left unchanged. This requires
[RDKit](https://www.rdkit.org/) (`pip install rdkit`); it is imported lazily, so RDKit
is only needed if you call this method.

### Removing heterogens — `removeHeterogens`

`removeHeterogens(keepWater=True, keepCoenzyme=False)` controls water and ligands
independently. Standard polymer residues and protein chain caps are always kept.

| Goal | Call |
| --- | --- |
| Drop water and ligands | `removeHeterogens(keepWater=False, keepCoenzyme=False)` |
| Drop ligands, keep water | `removeHeterogens(keepWater=True,  keepCoenzyme=False)` |
| Drop water, keep ligands | `removeHeterogens(keepWater=False, keepCoenzyme=True)` |
| Keep everything | `removeHeterogens(keepWater=True,  keepCoenzyme=True)` |

### Putting it together

A typical preparation that removes solvent but keeps and hydrogenates a known ligand:

```python
fixer = PDBFixer(pdbid='3PTB')
fixer.removeHeterogens(keepWater=False, keepCoenzyme=True)   # drop water, keep ligand
fixer.findMissingResidues(); fixer.missingResidues = {}
fixer.findMissingAtoms(); fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)                               # protein + caps
fixer.addHydrogensToLigands({'BEN': 'c1ccc(cc1)C(=N)N'})     # ligand
fixer.atomicOPT()                                            # energy-minimize (optional)
```

If you remove the ligands instead, you do not need a template or RDKit at all — just
`removeHeterogens(...)` followed by `addMissingHydrogens`.

## Installation

PDBFixer can be installed with conda or mamba.

```
conda install -c conda-forge pdbfixer
```

Alternatively you can install from source, as described in the manual.