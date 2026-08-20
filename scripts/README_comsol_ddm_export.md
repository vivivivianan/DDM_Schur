# Generic COMSOL DDM exporter

Start COMSOL 6.1 LiveLink/MATLAB, then invoke:

```matlab
export_comsol_ddm_case("my_model.mph", "my_partition.json", "output/my_case")
```

The JSON assigns actual COMSOL Domain IDs to subdomains and supplies each independent mesh size in metres. The exporter never saves or overwrites the input MPH. It creates temporary in-memory `ddm_mesh_###` mesh sequences, exports unchanged COMSOL MPHTXT files, validates their entity records, and writes the existing C++ solver input grammar.

Version 1 accepts constant isotropic `thermalconductivity`, `density`, and `heatcapacity`; Temperature, Convective Heat Flux, inward Heat Flux, Heat Source, Initial Values, and basic transient timing. Unsupported expressions, materials, or physics fail explicitly.
