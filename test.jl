include("lmpio.jl")
include("manyMolecules.jl")
include("packing.jl")

using .manyMolecules, .PackingMolecules, .LammpsIO

config = Vector{Any}(undef, 0)
#push!(config, (filename="molecules/benzene.data", molname="benzene", nmol=400))
push!(config, (filename="molecules/1-undecene.data", molname="1-undecene", nmol=100))
#push!(config, (filename="molecules/cis-2-butene.data", molname="butene", nmol=400))
#push!(config, (filename="molecules/cis-3-hexene.data", molname="cis-3-hexene", nmol=400))

atoms, topologies, coeffs = manyMolecules.add_molecule(config)

atoms, box = PackingMolecules.init_position(atoms)

name ="lmpsystem/1-undecene-100L1.data"
LammpsIO.write_lmpdat(name, atoms, topologies, coeffs, box)