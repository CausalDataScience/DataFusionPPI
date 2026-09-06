# DataFusionPPI

## Purpose

DataFusionPPI is a model-agnostic paper workspace for prediction-powered fusion of randomized-trial (RCT) and observational-study (OBS) data, targeting conditional average treatment effects (CATE) and average treatment effects (ATE) in the RCT population.

The canonical seed is `memo/PPI-Fusion.md`. It was copied verbatim from `research/papers/DataFusionPFN/study/PPI-Fusion.md` during project initialization.

## Project Boundary

- `DataFusionPPI` develops the model-agnostic PPI data-fusion formulation.
- `DataFusionPFN` is the sibling workspace for PFN-specific implementation and experiments.
- Claims, results, and implementation status are not transferred between the projects without explicit verification.

## Reconstruction Provenance

This workspace was freshly reconstructed on 2026-08-14 CDT. Preflight found no existing `research/papers/DataFusionPPI/` path and no Git history for that exact path. A similarly named nested draft surface under `DataFusionPFN` was not treated as prior history for this sibling project.

The manuscript was materialized from `research/templates/paper/latex_skeleton_itbound_struct_v2/manuscript/*.tpl`. Initialization changed only the paper title and template claim placeholders. The template text is scaffolding, not a verified novelty, theorem, or empirical-result statement for this project.

## Structure

- `memo/PPI-Fusion.md`: canonical seed
- `manuscript/`: materialized LaTeX scaffold
- `code/`: reserved for future approved implementation
- `materials/`: reserved for future approved source artifacts
- `STATUS.md`: current project state and verification boundary

## Current State

The project is waiting for substantive research specification and validation. No novelty claim, theorem, proof, estimator guarantee, simulation result, or empirical result has been verified or asserted during this structural setup.
