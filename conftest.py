"""Put ``src/`` on sys.path so ``import fluvtree`` resolves during test collection
without an editable install (src-layout). ``pip install -e .`` is the production
path; this keeps a bare ``pytest`` working in a fresh checkout."""

import os
import sys

_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
