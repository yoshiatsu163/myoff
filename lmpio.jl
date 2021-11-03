module LammpsIO
#export read_lmpdat, write_lmpdat

using Pipe, JuliaDB

mutable struct Box
    xdim::String
    ydim::String
    zdim::String
    triclinic::String
end

mutable struct Header
    comment::String
    n_attributes::Vector{Int64}
    n_types::Vector{Int64}
    box::Box
end

mutable struct Sections
    header::Vector{String}
    topologies::Vector{String}
    atoms::Vector{String}
    bonds::Vector{String}
    angles::Vector{String}
    diheds::Vector{String}
    imprs::Vector{String}

    Sections() = new()
end

function read_lmpdat(filename::String)
    fp=open(filename, "r")
    lines = readlines(fp)
    close(fp)

    sname = ["Atoms"     , "Bonds"     , "Angles"    , "Dihedrals" , "Impropers" ]
    field = Dict(
        sname[1] => :topologies,
        sname[2] => :atoms,
        sname[3] => :bonds,
        sname[4] => :angles,
        sname[5] => :diheds
    )

    sections = Sections()

    f = findfirst(line->occursin("Masses", line), lines)
    sections.header = @pipe filter(!isempty, lines[1:f-1]) |> map(l->replace(l, "\t"=>"    "), _)

    no_impr = false
    for str in sname
        s = f
        f = findfirst(line->occursin(str, line), lines)
        if f==nothing && str=="Impropers"
            no_impr = true
            f = length(lines)+1
            @pipe filter(!isempty, lines[s:f-1]) |> map(l->replace(l, "\t"=>"    "), _) |> setfield!(sections, field[str], _)
            break
        elseif f != nothing
            @pipe filter(!isempty, lines[s:f-1]) |> map(l->replace(l, "\t"=>"    "), _) |> setfield!(sections, field[str], _)
        else
            pritnln("Missing Atoms or Bonds or Angles or Dohedrals in the lmpdat.")
            exit()
        end
    end

    if no_impr
        sections.imprs = ["Impropers"]
    else
        sections.imprs = @pipe filter(!isempty, lines[f:end]) |> map(l->replace(l, "\t"=>"    "), _)
    end

    sections
end

function write_lmpdat(filename, atoms, topologies, coeffs, box)
    fp = open(filename, "w")

    #header
    write(fp, "# This data file is created by manyMolecule.jl\n\n")
    @pipe select(atoms, :atomid) |> maximum |> write(fp, "$_ atoms\n")
    @pipe filter(t->t[:category]=="Bonds"    , topologies) |> length(select(_, :tid)) |> write(fp, "$_ bonds\n")
    @pipe filter(t->t[:category]=="Angles"   , topologies) |> length(select(_, :tid)) |> write(fp, "$_ angles\n")
    @pipe filter(t->t[:category]=="Dihedrals", topologies) |> length(select(_, :tid)) |> write(fp, "$_ dihedrals\n")
    @pipe filter(t->t[:category]=="Impropers", topologies) |> length(select(_, :tid)) |> write(fp, "$_ impropers\n")
    write(fp, "\n")

    @pipe select(atoms, :atomtype) |> maximum |> write(fp, "$_ atom types\n")
    @pipe filter(t->t[:category]=="Bond Coeffs"    , coeffs) |> length(select(_, :cid)) |> write(fp, "$_ bond types\n")
    @pipe filter(t->t[:category]=="Angle Coeffs"   , coeffs) |> length(select(_, :cid)) |> write(fp, "$_ angle types\n")
    @pipe filter(t->t[:category]=="Dihedral Coeffs", coeffs) |> length(select(_, :cid)) |> write(fp, "$_ dihedral types\n")
    @pipe filter(t->t[:category]=="Improper Coeffs", coeffs) |> length(select(_, :cid)) |> write(fp, "$_ improper types\n")
    write(fp, "\n")
    
    write(fp, "0.0 $(box[1]) xlo xhi\n")
    write(fp, "0.0 $(box[2]) ylo yhi\n")
    write(fp, "0.0 $(box[3]) zlo zhi\n")
    write(fp, "0.0 0.0 0.0 xy xz yz\n")

    write(fp, "\nMasses\n\n")
    for molname in select(atoms, :molname) |> unique
        for row in filter(t->(t[:category]=="Masses" && t[:molname]==molname), coeffs)
            write(fp, "$(row[:cid])  $(row[:coeff])\n")
        end
    end
    write(fp, "\n")

    for cat in @pipe select(coeffs, :category) |> unique |> filter(!=("Masses"),_)
        write(fp, "$(cat)\n\n")
        for row in filter(t->t[:category]==cat, coeffs)
            write(fp, "$(row[:cid]) $(row[:style]) $(row[:coeff])\n")
        end
        write(fp, "\n")
    end

    write(fp, "Atoms\n\n")
    for row in atoms
        write(fp, "$(row[:atomid])  $(row[:molid])  $(row[:atomtype])    $(row[:q])     $(row[:x][1])     $(row[:x][2])     $(row[:x][3])\n")
    end
    write(fp, "\n")

    for cat in select(topologies, :category) |> unique
        write(fp, "$(cat)\n\n")
        for row in filter(t->t[:category]==cat, topologies)
            str = "$(row[:tid]) $(row[:cid])"
            for aid in row[:atomid]
                str *= "     $(aid)"
            end
            write(fp, str * "\n")
        end
        write(fp, "\n")
    end

    close(fp)
end

###
end

