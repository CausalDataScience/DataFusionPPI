"""Confirmatory CATE pipeline for DataFusionPPI.

This package is separate from the exploratory modules in ``code/`` on purpose.
The audit in ``materials/2026-09-12-cate-experiment-audit-and-confirmatory-plan.md``
found that the exploratory pipeline and the confirmatory requirements had drifted
apart, so nothing here imports the exploratory drivers, and every artifact it
writes lands under ``materials/confirm/`` with an explicit exploratory or
confirmatory label.
"""
