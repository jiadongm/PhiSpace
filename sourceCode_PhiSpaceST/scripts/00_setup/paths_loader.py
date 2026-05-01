"""PhiSpace ST path loader (Python).

Usage from any analysis script::

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "00_setup"))
    from paths_loader import paths

It walks up from the current working directory to find
``config/paths_template.yml``, reads ``paths_local.yml`` if present (else the
template, with a warning), and exposes a ``paths`` dict.

Convenience derivations::

    paths["data_root"]   = paths["phispace_data_root"] + "/data"
    paths["output_root"] = paths["phispace_data_root"] + "/output"
"""

from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path

try:
    import yaml
except ImportError as exc:
    raise ImportError(
        "paths_loader: the 'pyyaml' package is required. "
        "Install with `pip install pyyaml`."
    ) from exc


def _find_proj_root(start: str | os.PathLike | None = None) -> Path:
    p = Path(start or os.getcwd()).resolve()
    while not (p / "config" / "paths_template.yml").exists():
        if p.parent == p:
            raise FileNotFoundError(
                "paths_loader: could not find sourceCode_PhiSpaceST/config/"
                "paths_template.yml by walking up from '"
                f"{Path(start or os.getcwd()).resolve()}'."
            )
        p = p.parent
    return p


def _load() -> dict:
    proj = _find_proj_root()
    local_yml = proj / "config" / "paths_local.yml"
    template_yml = proj / "config" / "paths_template.yml"
    cfg_path = local_yml if local_yml.exists() else template_yml
    if cfg_path == template_yml:
        warnings.warn(
            "[paths_loader] Using paths_template.yml (placeholders). "
            "Copy it to paths_local.yml and edit before running real analyses.",
            stacklevel=2,
        )
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f) or {}
    cfg["project_root"] = str(proj)
    pdr = cfg.get("phispace_data_root")
    if pdr:
        cfg.setdefault("data_root", str(Path(pdr) / "data"))
        cfg.setdefault("output_root", str(Path(pdr) / "output"))
    return cfg


paths = _load()
