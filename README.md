# WignerSeitzDefectAnalysis

`WignerSeitzDefectAnalysis` identifies vacancies, interstitials and antisites by comparing
a current atomic configuration against a reference lattice using Wigner-Seitz cells.

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- WignerSeitzDefectAnalysis
```

## Build

```bash
conan install . --build=missing -s compiler.cppstd=17 -o boost/*:without_stacktrace=True
cmake --preset conan-release
cmake --build --preset conan-release
```

## CLI

Usage:

```bash
wigner-seitz-defect-analysis <lammps_file> [output_base] --reference <ref_file> [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file (current configuration). | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--reference <path>` | Yes | Reference LAMMPS dump file defining the lattice sites. | |
| `--affineMapping <mode>` | No | Affine mapping mode: `off`, `toReference`, `toCurrent`. | `off` |
| `--eliminateCellDeformation` | No | Eliminate cell deformation before site assignment (currently a no-op in the engine). | `false` |
| `--minimumImageConvention <bool>` | No | Use minimum image convention for site assignment. | `true` |
| `--perTypeOccupancies` | No | Emit occupancy counts split by particle type. | `false` |
| `--threads <int>` | No | Maximum worker threads. | auto |
| `--help` | No | Print CLI help. | |

> The current and reference configurations may (and usually do) have different atom counts -
> that mismatch is how vacancies and interstitials are detected.

## Outputs

Two Parquet files are produced under `{output_base}_*.parquet`:

- `{output_base}_wigner_seitz.parquet` - per-site / summary defect listing with
  `main_listing` (vacancy / interstitial / antisite / occupied counts, total sites)
  and a `per-site-properties` array.
- `{output_base}_atoms.parquet` - per-role records (the canonical per-atom schema,
  D-009). The rows are the reference sites plus extra interstitial occupants (a
  synthetic frame), so vacancies - which have no current atom - are represented.
  The `bucket` column holds the role (`Vacancy`, `Occupied`, `Interstitial`,
  `Antisite`) and drives the `AtomisticExporter` GLB grouping. Per-atom columns:
  `occupancy`, `site_index`, and `occupancy_per_type` (when `--perTypeOccupancies`).
