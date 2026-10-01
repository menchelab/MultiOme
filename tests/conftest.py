from __future__ import annotations

import numpy as np
import pytest

from multiome_core.schema import GeneGroup, Layer, Multiplex


def planted_multiplex(seed: int = 0, n: int = 300, module: int = 30):
    """Three random layers; genes G000..G029 form a dense module in layers A and B only."""
    rng = np.random.default_rng(seed)
    genes = [f"G{i:03d}" for i in range(n)]
    mod = genes[:module]

    def layer(lid, p_bg, planted, drop):
        src, dst = [], []
        keep = [g for g in genes if g not in drop]
        for _ in range(int(p_bg * len(keep))):
            a, b = rng.choice(len(keep), 2, replace=False)
            src.append(keep[a])
            dst.append(keep[b])
        if planted:
            for i in range(module):
                for j in range(i + 1, module):
                    if rng.random() < 0.3:
                        src.append(mod[i])
                        dst.append(mod[j])
        return Layer.from_edges(lid, src, dst, tags=(lid.lower(),))

    mpx = Multiplex(name="planted")
    mpx.add(layer("A", 3, True, set()))
    mpx.add(layer("B", 3, True, set(genes[250:])))  # B lacks 50 genes
    mpx.add(layer("C", 3, False, set()))
    group = GeneGroup("module", mod, label="Planted module")
    others = [GeneGroup(f"rand{k}", list(rng.choice(genes[module:], 25, replace=False)))
              for k in range(2)]
    return mpx, [group, *others]


@pytest.fixture(scope="session")
def planted():
    return planted_multiplex()
