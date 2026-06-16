# Wigner-Seitz Defect Analysis

Identifies vacancies, interstitials and antisites by comparing the current configuration against a reference lattice using Wigner-Seitz cells.

## Install

```bash
vpm install @voltlabs/wigner-seitz-defect-analysis
```

## CLI

```bash
wigner-seitz-defect-analysis <input_dump> [output_base] --reference <ref_file> [options]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `<input_dump>` | yes | — | Input LAMMPS dump (current configuration). |
| `[output_base]` | no | derived from input | Base path for output files. |
| `--reference <path>` | yes | — | Reference LAMMPS dump defining the lattice sites. |
| `--affineMapping <mode>` | no | `off` | Affine mapping mode: `off`, `toReference`, `toCurrent`. |
| `--eliminateCellDeformation` | no | `false` | Eliminate cell deformation before site assignment. |
| `--minimumImageConvention <bool>` | no | `true` | Use minimum image convention for site assignment. |
| `--threads <int>` | no | auto | Maximum worker threads. |
| `--help` | no | — | Print CLI help. |

## Exports

| Output file | Exposure | Exporter → artifact |
|---|---|---|
| `{output_base}_wigner_seitz.parquet` | Wigner-Seitz Sites | — |
| `{output_base}_atoms.parquet` | Wigner-Seitz Model | AtomisticExporter → glb |

---

Full input contract and examples: https://docs.voltcloud.dev/docs/plugins
