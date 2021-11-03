include("lmpio.jl")
using .LammpsIO

module ManyMolecules

using StaticArrays, Statistics
import Base.+

for str in ["Atom", "Bond", "Dihed", "Impr"]
    for suf in ["id", "Type"]
        Meta.parse("struct $(str)$(suf)
            toi::Int64
            $(str)$(suf)(toi::Int64) = new(toi)
            $(str)$(suf)(str::String) = new(parse(Int64, str))
        end"
        ) |> eval

        Meta.parse("function +(n::$(str)$(suf), m::$(str)$(suf))
            $(str)$(suf)( n.toi + m.toi )
        end"
        ) |> eval
    end
    #Meta.parse("struct $(str)id toi::Int64 end") |> eval
    #Meta.parse("struct $(str)Type toi::Int64 end") |> eval
    Meta.parse("mutable struct $(str)coeff  \nid::$(str)Type  \nrest::String  \nend") |> eval
end
struct MoleculeId toi::Int64 end

mutable struct Paircoeff
    id::AtomType
    rest::String
end

struct Charge
    tof::Float64
end

struct Position
    unwrap::SVector{3, Float64}
end

mutable struct Atom
    id::Atomid
    type::AtomType
    charge::Charge
    x::Position
end

mutable struct Bond
    id::Bondid
    type::BondType
    origin::Atomid
    target::Atomid
end

mutable struct Angle
    id::Angleid
    type::AngleType
    head::Atomid
    center::Atomid
    tail::Atomid
end

mutable struct Dihedral
    id::Dihedid
    type::DihedType
    head::Atomid
    c1::Atomid
    c2::Atomid
    rear::Atomid
end

mutable struct Impr
    id::Imprid
    type::ImprType
    a1::Atomid
    a2::Atomid
    a3::Atomid
    a4::Atomid
end

mutable struct BoundingCircle
    center::SVector{3, Float64}
    radius::Float64
end
function BoundingCircle(mol::Molecule)
    N = length(mol.atoms)
    centroid = [mean.(mol.atoms[i]) ./N for i in 1:3]
    radius = 0.0
    for vec in mol.atoms
        rsq = (vec .- centroid) .^2 |> sum
        ifelse(radius < rsq, radius = rsq, continue)
    end

    BoundingCircle(centroid, radius)
end

mutable struct Molecule
    id::MolculeId
    atoms::Vector{Atom}
    bonds::Vector{Bond}
    angles::Vector{Angle}
    diheds::Vector{Dihedral}
    imprs::Vector{Impr}
    bcircle::BoundingCircle

    Molecule() = new()
end
function Molecule(sections, molid)
    molecule = Molecule()
    molecule.id = molid

    for arr in split.(filter(!=(""), sections.atom)) 
        push!(molecule.atoms, 
            Atom(
                parse(Int64, arr[1]) |> Atomid, 
                parse(Int64, arr[2]) |> AtomType, 
                parse(Float64, arr[4]) |> Charge, 
                parse.(Float64, arr[5:end]) |> Position, 
            )
        )
    end

    for arr in split.(filter(!=(""), sections.bond)) 
        push!(molecule.bonds, 
            Bond(
                parse(Int64, arr[1]) |> Bondid, 
                parse(Int64, arr[2]) |> BondType, 
                parse(Int64, arr[3]) |> Atomid, 
                parse(Int64, arr[4]) |> Atomid, 
            )
        )
    end

    for arr in split.(filter(!=(""), sections.angle)) 
        push!(molecule.angles, 
            Angle(
                parse(Int64, arr[1]) |> Angleid, 
                parse(Int64, arr[2]) |> AngleType, 
                parse(Int64, arr[3]) |> Atomid, 
                parse(Int64, arr[4]) |> Atomid, 
                parse(Int64, arr[5]) |> Atomid, 
            )
        )
    end

    for arr in split.(filter(!=(""), sections.diheds)) 
        push!(molecule.diheds, 
            Dihedral(
                parse(Int64, arr[1]) |> Dihedid, 
                parse(Int64, arr[2]) |> DihedType, 
                parse(Int64, arr[3]) |> Atomid, 
                parse(Int64, arr[4]) |> Atomid, 
                parse(Int64, arr[5]) |> Atomid, 
                parse(Int64, arr[6]) |> Atomid, 
            )
        )
    end
       
    for arr in split.(filter(!=(""), sections.imprs)) 
        push!(molecule.imprs, 
            Impr(
                parse(Int64, arr[1]) |> Imprid, 
                parse(Int64, arr[2]) |> ImprType, 
                parse(Int64, arr[3]) |> Atomid, 
                parse(Int64, arr[4]) |> Atomid, 
                parse(Int64, arr[5]) |> Atomid, 
                parse(Int64, arr[6]) |> Atomid, 
            )
        )
    end

    molecule.bcircle = BoundingCircle(molecule)
    molecule
end

mutable struct Substance
    molecules::Vector{Molecule}
    paircoeffs::Vector{Paircoeff}
    bondcoeffs::Vector{Bondcoeff}
    anglecoeffs::Vector{Anglecoeff}
    dihedcoeffs::Vector{Dihedcoeff}
    imprcoeffs::Vector{Imprcoeff}

    Substance() = new()
end
function Substance(sections::Sections)
    keys = ["Pair", "Bond", "Angle", "Dihedral", "Improper", "Atoms"]
    n = [findfirst(line->occursin(key, line), sections.header) for key in keys]
    substance = Substance()

    #arr = @pipe sections.header[n[1]:n[2]] |> filter(l->occursin(r"^[1-9]", l), _) |> split.(_)
    #substance.paircoeffs = Vector{PairCoeff}(undef, length(arr))
    #for j in 1:length(arr)
    #    substance.paircoeffs[j].id = Atomid(arr[j][1])
    #    substance.paircoeffs[j].rest = join(arr[j][2:end], " ")
    #end

    names = ["Atom", "Bond", "Angle", "Dihed", "Impr"]
    dict = Dict("Atom"=>"pair", "Bond"=>"bond", "Angle"=>"angle", "Dihed"=>"dihed", "Impr"=>"impr")
    for i in 1:length(keys)-1
        typename = ""
        ifelse(names[i]=="Atom", typename=="Pair", typename=names[i])
        idtype = names[i]
        coeffname = dict[idtype]
        "arr = @pipe sections.header[n[i]:n[i+1]] |> filter(l->occursin(r\"^[1-9]\", l), _) |> split.(_)
        substance.$(coeffname)coeffs = Vector{$(typeame)Coeff}(undef, length(arr))
        for j in 1:length(arr)
            substance.$(coeffname)coeffs[j].id = $(idtype)id(arr[j][1])
            substance.$(coeffname)coeffs[j].rest = join(arr[j][2:end], " ")
        end" |> eval
    end

    substance
end

function minimalSubstance(sections)
    substance = Substance(sections)
    push!(
        substance.molecule,
        Molecule(sections, MoleculeId(1))
    )

    substance
end

function create_molecule_sets(config::Vector{NamedTuple{(:substance, :nmol), Tuple{Substance, Int64}}})
    @assert mapreduce(x -> x.nmol>0, &, config) "Number of any molecule must be >0"

    substances = Vector{Substance}(undef, 0)
    for i in 1:mapreduce(x->x.nmol, *, config)
        push!(substances, substance_with_unique_id!(config, substances))
    end
    uniformize_atomid!(substances)

    substances
end

function substance_with_unique_id!(config, substances)
    next = copy(config[end].substance)
    nmol = config[end].nmol
    nexyt.molecules = [ next.molecules for i in 1:nmol ]
    for i in 1:nmol
        next.molecules[i].id = MoleculeId(i)
    end
    if length(substances) == 0
        pop!(config)
        return next
    end
    
    prev = substances[end]
    @assert fieldnames(typeof(next)) == (:molecules, :paircoeffs, :bondcoeffs, :anglecoeffs, :dihedcoeffs, :imprcoeffs)
    f2f = Dict( :paircoeffs=>:atoms, :bondcoeffs=>:bonds, :anglecoeffs=>:angles, :dihedcoeffs=>:diheds, :imprcoeffs=>:imprs )
    for sym in fieldnames(typeof(next))
        new_coeffid = Dict()
        next_field = getfield(next, sym)
        prev_field = getfield(prev, sym)
        for i in 1:length(next_field)
            offset_id = next_field[i].id + prev_field[end].id
            push!(new_coeffid, next_field[i].id => offset_id)
            next_field[i].id = offset_id
        end
        sym == :molecules && continue

        for i in 1:length(next.molecules)
            atom_likes = getfield(next.molecules[i], f2f[sym])
            for j in 1:length(atom_likes)
                atom_likes[j].type = newid[atom_likes[j].type]
            end
        end
    end
    pop!(config)
    next
end

function uniformize_atomid!(substances)
    natom = mapreduce(x -> length(x.molecules)*length(x.molecules[1].atoms), +, substances)
    cnt = Atomid(1)
    for i in 1:length(substances), j in 1:length(substances[i].molecules)
        new_atomid = Dict()
        for sym in fieldnames(Molecule) |> x->filter(!=(:id), x) |> x->filter(!=(:bcircle), x)
            for k in 1:length(substances[i].molecules[1].atoms)
                push!(newid, substances[i].molecules[j].atoms[k].id => cnt)
                substances[i].molecules[j].atoms[k].id = cnt
                cnt += Atomid(1)
            end
            for k in 1:length(substances[i].molecules[1].bonds)
                substances[i].molecules[j].bonds[k]. = 
            end
    end

    for i in length(substances[i].molecules[1].bonds)
        substances[i].molecules[1].bonds
end

struct Box
    x::SVector{3, Float64}
    y::SVector{3, Float64}
    z::SVector{3, Float64}
end

mutable struct MDSystem
    substances::Vector{Substance}
    atom2mass::Dict{Atomid, Foat64}
    box::Box
end
function MDSystem(substances::Vector{Substance}, tolerance=3.5)
    #attribute id の重複を排除
    #atom idの重複を排除
    rm_id_duplication!(substances)
        
    sunstance.molecule.atom
    #分子が干渉しないよう平行移動
    #bonの設定
end





##########################################################
end


