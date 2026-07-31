"""
Common physics modules: optional, opt-in physics reusable across solution types.

Nothing here is required to run the core (network + diffusion solve + basic
transport); a model *includes* a common module when it wants that physics, and
different solvers -- and, eventually, external models like GRLP -- can share the
same implementation from here. Each module is off by default.

First member: :mod:`fluvtree.common.gravel_attrition` -- Sternberg downstream
fining of the transported gravel load (GRLP's ``update_gravel_loss``), which
applies to any river carrying gravel, gravel-bedded or bedrock-floored.
"""
