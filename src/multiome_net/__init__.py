"""multiome_net — the network-generation unit (phase 2+).

Fetches data (ToolUniverse, STRING, DepMap, GTEx, ...) and builds network layers,
emitting standardized .mpx bundles consumed by multiome_algo. Layers are pluggable:
co-expression (correlation or regulatory inference), PPI, ontology similarity, etc.

Phase 1 uses multiome_core.legacy instead; this package is scaffolding for now.
"""
