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

### Which method hydrogenates what

| Method | Hydrogenates | Needs network / CCD |
| --- | --- | --- |
| `addMissingHydrogens(pH=7.0)` | protein residues **and chain caps** (`ACE`, `NME`, `NH2`, `NMA`, `FOR`); ligands left untouched | no (for standard residues) |
| `addHydrogensToLigands(templates)` | ligands/heterogens you supply a SMILES (or RDKit `Mol`) for | **no — fully offline** |

The two are independent. Call only `addMissingHydrogens` to protonate the protein and
leave ligands bare; add `addHydrogensToLigands` when you also want the ligands
hydrogenated.

### Hydrogenating the protein — `addMissingHydrogens`

```python
from pdbfixer import PDBFixer

fixer = PDBFixer(filename='protein.pdb')
fixer.findMissingResidues(); fixer.missingResidues = {}   # skip loop building
fixer.findMissingAtoms(); fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)                            # protein + chain caps
```

Standard amino acids and the chain caps `ACE`/`NME`/`NH2`/`NMA`/`FOR` are hydrogenated;
ligands, ions, and other heterogens are carried over unchanged (no hydrogens added).

### Hydrogenating ligands with RDKit — `addHydrogensToLigands`

The hydrogen topology of an arbitrary small molecule **cannot** be recovered reliably
from heavy-atom coordinates alone (bond orders are ambiguous). So instead of guessing,
you supply the chemistry you already know — a SMILES string (or RDKit `Mol`) per ligand
residue name. RDKit transfers its bond orders onto the matching heavy atoms and adds the
hydrogens. It runs **fully offline**: no Chemical Component Dictionary, no network.

```python
fixer = PDBFixer(filename='complex.pdb')
fixer.addMissingHydrogens(7.0)                                    # protein + caps
added = fixer.addHydrogensToLigands({'LIG': 'CC(=O)Oc1ccccc1C(=O)O'})
# added == {'LIG': 8}  -> 8 hydrogens were added to residue LIG
```

Multiple ligands at once — only the residues you provide a SMILES for are touched:

```python
fixer.addHydrogensToLigands({
    'BEN': 'c1ccc(cc1)C(=N)N',     # benzamidine
    'LIG': 'CC(=O)Oc1ccccc1C(=O)O' # aspirin
})
```

Behaviour and guarantees:

- **Heavy-atom positions are never changed.** Only the new hydrogens get coordinates;
  every existing heavy atom keeps its exact input position.
- The template value may be a **SMILES string** or an **RDKit `Mol`**. Hydrogens in the
  template are ignored.
- The template must contain **exactly the same heavy atoms** (same element graph) as the
  residue in the model. If the count/graph does not match, that ligand is **left
  unchanged** and a warning is printed — the structure is never corrupted.
- A ligand with **no** entry in `templates` is left unchanged. This is the control:
  only ligands you name (and provide a SMILES for) are hydrogenated.
- Returns a dict mapping residue name -> number of hydrogens added.
- Requires [RDKit](https://www.rdkit.org/) (`pip install rdkit`); it is imported lazily,
  so RDKit is only needed if you actually call this method.

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

Remove water, keep the ligand, hydrogenate protein **and** ligand, then minimize:

```python
fixer = PDBFixer(pdbid='3PTB')
fixer.removeHeterogens(keepWater=False, keepCoenzyme=True)   # drop water, keep ligand
fixer.findMissingResidues(); fixer.missingResidues = {}
fixer.findMissingAtoms(); fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)                               # protein + caps
fixer.addHydrogensToLigands({'BEN': 'c1ccc(cc1)C(=N)N'})     # ligand (heavy atoms unmoved)
fixer.atomicOPT()                                            # energy-minimize (optional)
```

Remove water **and** ligands, then hydrogenate the protein only (no SMILES / RDKit
needed, because nothing but protein remains):

```python
fixer = PDBFixer(pdbid='3PTB')
fixer.removeHeterogens(keepWater=False, keepCoenzyme=False)  # drop water + ligands
fixer.findMissingResidues(); fixer.missingResidues = {}
fixer.findMissingAtoms(); fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)                              # protein + caps
```

Write the result out (PDB or PDBx/mmCIF):

```python
from openmm.app import PDBFile, PDBxFile
PDBFile.writeFile(fixer.topology, fixer.positions, open('out.pdb', 'w'))
# or:  PDBxFile.writeFile(fixer.topology, fixer.positions, open('out.cif', 'w'))
```

## Installation

PDBFixer can be installed with conda or mamba.

```
conda install -c conda-forge pdbfixer
```

Alternatively you can install from source, as described in the manual.