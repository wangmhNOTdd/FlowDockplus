# Project Guidelines

## Code Style
- Python style/quality is enforced through `pre-commit` hooks. Run `make format` before提交（见 `.pre-commit-config.yaml`, `Makefile`）。
- Test config comes from `pyproject.toml` (`pytest` markers, doctest, coverage exclusions).
- Follow existing entrypoint style: Hydra + Lightning + `rootutils.setup_root(...)` (`flowdock/train.py`, `flowdock/eval.py`, `flowdock/sample.py`).
- Keep edits small and local; avoid refactoring unrelated modules in this research codebase.

## Architecture
- Main runtime entrypoints:
  - Training: `flowdock/train.py`
  - Evaluation: `flowdock/eval.py`
  - Sampling/inference: `flowdock/sample.py`
- Core model logic lives in `flowdock/models/flowdock_fm_module.py`.
- Data orchestration is in `flowdock/data/combined_datamodule.py` with dataset components under `flowdock/data/components/`.
- Config graph is Hydra-based under `configs/` with task roots: `train.yaml`, `eval.yaml`, `sample.yaml`.

## Build and Test
- Environment (recommended):
  - `mamba env create -f environments/flowdock_environment.yaml`
  - `conda activate FlowDock`
  - `pip install -e .`
- Quality checks:
  - `make format`
- Tests:
  - `make test` (exclude `slow`)
  - `make test-full`
- Run tasks:
  - `python flowdock/train.py trainer=gpu`
  - `python flowdock/eval.py ckpt_path=<path> trainer=gpu`
  - `python flowdock/sample.py ckpt_path=<path> out_path=<path> ...`

## Project Conventions
- Use Hydra overrides instead of hardcoding values. Typical pattern: `python flowdock/train.py experiment=flowdock_fm trainer.max_epochs=20`.
- Paths are resolved via project root (`.project-root`) and `configs/paths/default.yaml`; prefer config-driven paths.
- Sampling supports either direct inputs (`input_receptor`, `input_ligand`) or batched CSV (`csv_path`) with columns: `id,input_receptor,input_ligand,input_template` (`flowdock/sample.py`, `configs/sample.yaml`).
- For strategy/environment selection, rely on config groups (`configs/strategy/*`, `configs/environment/*`) rather than inline trainer rewrites.

## Sampling Pipeline: Template and Prior Logic
The `predict_step` in `flowdock/models/flowdock_fm_module.py` follows this priority:

1. **`input_template` provided and file exists** → uses the PDB directly as `apo_rec_path` (the VD-ODE starting structure). ESMFold weights are **not** loaded (skipped for efficiency).
2. **`input_template=null` + `prior_type=esmfold`** → calls `esmfold_model.infer_pdb()` at runtime to predict apo structure from sequence. ESMFold weights **are** loaded.
3. `input_receptor` as a sequence string → a zero-coordinate dummy PDB is created (`create_full_pdb_with_zero_coordinates`) solely for sequence/ESM2 embedding extraction; actual 3D coordinates come from `apo_rec_path`.

**Template file naming conventions in `data/posebusters_benchmark_set/`:**
- `posebusters_benchmark_holo_aligned_esmfold_predicted_structures/` — ESMFold-predicted apo structures **pre-aligned to the holo crystal coordinate frame** (via `esmfold_apo_to_holo_alignment.py`). Use for benchmark RMSD evaluation.
- `posebusters_benchmark_esmfold_predicted_structures/` — raw ESMFold predictions, not aligned. Use for blind prediction.

**Benchmark vs blind prediction:**
- Reproducing paper results: use `holo_aligned` templates + `use_template=true`.
- True blind prediction (no crystal structure known): set `input_template=null`; ESMFold runs on-the-fly.

When providing `input_template`, the sequence extracted from the template PDB must exactly match `input_receptor`; a mismatch causes the sample to be skipped with an error log.

## Integration Points
- Core dependencies: PyTorch Lightning, Hydra/OmegaConf, RDKit, ESM/ESMFold-related processing.
- Logging/callback instantiation is centralized in `flowdock/utils/instantiators.py` and utility wrappers in `flowdock/utils/utils.py`.
- Optional cluster/HPC settings are under `configs/environment/slurm.yaml`.

## Security
- Do not hardcode secrets/tokens in configs or scripts.
- Treat large external downloads/checkpoints as untrusted inputs; keep extraction paths explicit (`README.md`, `checkpoints/`).
- GPU-memory-sensitive stages (sampling, ESMFold processing) should prefer config knobs (`chunk_size`, `esmfold_chunk_size`) over code changes.

## Current Research Direction (PhysFlowDock idea)
- Active proposal is documented in `idea.md`: inject chance-constrained projection into FlowDock VD-ODE sampling.
- Prefer minimal-invasive changes in inference path first (`flowdock/sample.py` + model sampling loop), not training pipeline rewrites.
- Keep this separation explicit in PRs: (1) constraint/scheduler/projection utilities, (2) sampler loop integration, (3) ablation toggles via config.
- Preserve baseline reproducibility: any new projection behavior must be gated by config flags and default to current FlowDock behavior.
