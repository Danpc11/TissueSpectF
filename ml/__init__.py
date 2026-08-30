"""TissueSpect-AE: a learned layer over the TissueSpectF statistical pipeline.

It does not replace anything. The statistical results — permutation p-values,
consensus spectra, meta-consensus — remain the evidence. This module learns a
representation of the same spectra and proposes components; a component only
this module finds is reported as an AI candidate, never as a finding.
"""

__version__ = "0.1.0"
