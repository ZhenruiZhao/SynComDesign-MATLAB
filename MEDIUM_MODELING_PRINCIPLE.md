# Medium Modeling Principle

The fixed MATLAB community model uses three explicit layers:

1. External medium exchange: `R_EX_*_u` reactions connect the outside world to the shared extracellular pool.
2. Shared extracellular pool: shared metabolites such as `ac[u]` or `no3[u]` are mass-balanced community metabolites.
3. Strain interface and internal metabolism: strain-prefixed exchange/interface and transport reactions connect each strain to the shared pool and its own compartments.

The medium file controls only external medium exchange. If a compound is absent from medium, its external shared exchange lower bound is closed to 0, so the outside world cannot supply it.

Cross-feeding is still allowed: if strain A secretes a compound into the shared pool, strain B can uptake it through its strain interface. Shared metabolite mass balance prevents strain B from consuming a compound that neither the medium nor another strain supplies.

Internal transport reactions are not controlled directly by medium. The previous community behavior was mixed because exchange detection was reused after prefixing and could close strain interfaces or transport-like reactions. The fixed mode avoids that by classifying community reactions and applying medium only to `external_medium_exchange`.