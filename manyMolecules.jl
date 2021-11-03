module manyMolecules

using ..LammpsIO
using JuliaDB, Pipe, StaticArrays, Setfield

struct Config
    filename::String
    molname::String
    nmol::Int64
end

mutable struct Offset
    molid::Int64
    atomid::Int64
    atomtype::Int64
    tid::Dict{String, Int64}
    cid::Dict{String, Int64}
end

mutable struct Topdict
    pair::Dict{Int64, Int64}
    bond::Dict{Int64, Int64}
    angle::Dict{Int64, Int64}
    dihed::Dict{Int64, Int64}
    impr::Dict{Int64, Int64}

    Topdict() = new(Dict(), Dict(), Dict(), Dict(), Dict())
end
    

const anames = [:molname, :molid, :atomid, :atomtype, :q, :x]
const tnames = [:molname, :molid, :category, :tid, :cid, :atomid]
const cnames = [:molname, :style, :category, :cid, :coeff]
const empty_atoms = table( 
            String[], Int64[], Int64[], Int64[], Float64[], Vector{Vector{Float64}}(undef, 0),
            names = anames
        )
const empty_topologies = table( 
                String[], Int64[], String[], Int64[], Int64[], Vector{Vector{Int64}}(undef, 0),
                names = tnames
            )
const empty_coeffs = table( 
            String[], String[], String[], Int64[], String[],
            names = cnames
        )

mutable struct Template
        atoms::Any
        topologies::Any
        coeffs::Any
        
        Template() = new(empty_atoms, empty_topologies, empty_coeffs)
        Template(atoms, topologies, coeffs) = new(atoms, topologies, coeffs)
end
function Base.copy(s::Template)
    p = select(s.atoms, :x)
    newa = table(
            select(s.atoms, :molname ) |> copy,
            select(s.atoms, :molid   ) |> copy,
            select(s.atoms, :atomid  ) |> copy,
            select(s.atoms, :atomtype) |> copy,
            select(s.atoms, :q       ) |> copy,
            copy.(select(s.atoms, :x)),
            #[ [p[i][1], p[i][2], p[i][3]] for i in 1:length(p) ], 
            names = anames
        )
    a = select(s.topologies, :atomid)
    newt = table(
            select(s.topologies, :molname ) |> copy,
            select(s.topologies, :molid   ) |> copy,
            select(s.topologies, :category  ) |> copy,
            select(s.topologies, :tid) |> copy,
            select(s.topologies, :cid) |> copy,
            copy.(select(s.topologies, :atomid)),
            #[ [a[i][j] for j in 1:length(a[i])] for i in 1:length(a) ],
            names = tnames
        )
    Template(newa, newt, s.coeffs)
end

const cat2cat = Dict(
    "Atoms"     => "Pair Coeffs",
    "Bonds"     => "Bond Coeffs",
    "Angles"    => "Angle Coeffs",
    "Dihedrals" => "Dihedral Coeffs",
    "Impropers" => "Improper Coeffs"
)

function add_molecule(config)
    if !is_molname_uniform(config)
        println("molname must be uniform.")
        return nothing
    end
    atoms = empty_atoms
    topologies = empty_topologies
    coeffs = empty_coeffs

    for conf in config
        sections = LammpsIO.read_lmpdat(conf.filename)
        molname = conf.molname
        nmol = conf.nmol
        # テーブルに新たな分子情報を追加，idを一意化する
        atoms, topologies, coeffs = push_tables!(atoms, topologies, coeffs, sections, molname, nmol)
    end
    #println(atoms)
    #println(topologies)
    #println(coeffs)

    return atoms, topologies, coeffs
end

function is_molname_uniform(config)
    molnames = [config[i].molname for i in 1:length(config)]
    unique(molnames) == molnames
end

function push_tables!(atoms, topologies, coeffs, sections, molname, nmol)
    #既存の一意化，ソート済みを確認
    detection = detect_non_uniform(atoms, topologies, coeffs)
    if detection != nothing
        println("Some duplications are detected at $(detection). There would be some bugs in push_table!().")
        println("A dump file \"tabledump\" is generated.")
        return nothing
    end

    # 最大値をoffsetとして取得
    offset = Offset(
        0, 0, 0, 
        Dict("Bonds"=>0, "Angles"=>0, "Dihedrals"=>0, "Impropers"=>0),
        Dict("Masses"=>0, "Pair Coeffs"=>0, "Bond Coeffs"=>0, "Angle Coeffs"=>0, "Dihedral Coeffs"=>0, "Improper Coeffs"=>0)
    )
    if length(atoms)!=0 && length(topologies)!=0 && length(coeffs) != 0
        offset.molid  = select(atoms, :molid)  |> maximum
        offset.atomid = select(atoms, :atomid) |> maximum
        offset.atomtype = select(atoms, :atomtype) |> maximum
        offset.tid    = Dict(cat => myfilter(topologies, (:category, ==(cat))) |> t->maximum(select(t, :tid)) for cat in select(topologies, :category) |> unique)
        offset.cid    = Dict(cat => myfilter(coeffs    , (:category, ==(cat))) |> t->maximum(select(t, :cid)) for cat in select(coeffs    , :category) |> unique)
    elseif length(atoms)==length(topologies)==length(coeffs) == 0
        
    else
        println("fatal error")
        exit()
    end

    #新規用の単分子テンプレートを作成
    template = create_template(molname, sections)
    #テンプレートを複製，一意化，テスト
    new = addend(template, offset, nmol)
    #データ結合
    atoms      = merge(atoms, new.atoms)
    topologies = merge(topologies, new.topologies)
    coeffs     = merge(coeffs, new.coeffs)

    atoms, topologies, coeffs
end

function detect_non_uniform(atoms, topologies, coeffs)
    # topologiesのatom id を unique |> sort してatomsと一致するか確認, 所属するmolid, molnameとの照合も行う

    if length(atoms)==length(topologies)==length(coeffs) == 0
        return nothing
    end

    #molid check
    ua = sort(select(atoms, :molid))
    ut = sort(select(topologies, :molid))
    @assert minimum(ua) == minimum(ut) == 1
    @assert ua == select(atoms, :molid) && ut == select(topologies, :molid) && unique(ua) == unique(ut)

    #atomtype check
    ua = select(atoms,  :atomtype)
    uc = myfilter(coeffs, (:category, ==("Pair Coeffs")) ) |> t->map(t->t[:cid], t)
    @assert minimum(ua) == minimum(uc) == 1
    @assert unique(ua) == sort(ua) |> unique "$(ua)"
    @assert uc==unique(uc)==sort(uc)
    @assert unique(ua)==uc

    # atomid check
    ids = select(atoms, :atomid)
    @assert ids == sort(ids)
    @assert unique(ids) |> sort == unique(ids)
    @assert minimum(ids) == 1

    #@pipe filter(t->t[:category]=="Bonds", topologies) |> select(_, :tid) |> println

    #tid check
    for cat in select(topologies, :category) |> unique
        tids = @pipe filter(t->t[:category]==cat, topologies) |> select(_, :tid)
        @assert unique(tids) == unique(tids) |> sort "$(cat): $(unique(tids))"
        @assert minimum(tids) == 1
    end

    # cid check
    categories = select(coeffs, :category) |> unique
    for catname in categories
        # 辞書が必要 "Bond Coeffs" => "Bonds"
        cid_coeffs = myfilter(coeffs,     (:category, ==(catname)) ) |> t->map(t->t[:cid], t)
        #cid_topols = myfilter(topologies, (:category, ==(catname)) ) |> t->map(t->t[:cid], t)
        @assert unique(cid_coeffs) |> sort == sort(cid_coeffs)
        #@assert unique(cid_topols) |> sort == sort(cid_topols)
        #@assert cid_topols[1] == cid_coeffs[1] == 1
    end

    nothing
end


function create_template(molname, sections)
    atoms = table( String[], Int64[], Int64[], Int64[], Float64[], Vector{Vector{Float64}}(undef, 0), names=anames )
    for i in 2:length(sections.atoms)
        arr = split(sections.atoms[i])
        tmp = table( [molname], [1], [parse(Int64,arr[1])], [parse(Int64,arr[3])], [parse(Float64, arr[4])], [parse.(Float64, arr[5:7])], names=anames )
        #append!( rows(atoms), rows(tmp) )
        atoms = merge(atoms, tmp)
    end

    topologies = table( String[], Int64[], String[], Int64[], Int64[], Vector{Vector{Int64}}(undef, 0), names=tnames )
    for sym in [:bonds, :angles, :diheds, :imprs]
        #getfield(sections, :diheds) |> println
        sec = getfield(sections, sym)
        for i in 2:length(sec)
            arr = split(sec[i])
            tmp = table( [molname], [1], [sec[1]], [parse(Int64,arr[1])], [parse(Int64,arr[2])], [parse.(Int64, arr[3:end])], names=tnames )
            #append!( rows(topologies), rows(tmp) )
            topologies = merge(topologies, tmp)
        end
    end

    # 名前逆に注意
    category = ""
    coeffs = table( String[], String[], String[], Int64[], String[], names=cnames )
    for i in 1:length(sections.topologies)
        if occursin(r"^[A-Z]", sections.topologies[i])
            category = sections.topologies[i]
            continue
        end
        arr = split(sections.topologies[i])
        if tryparse(Float64, arr[2]) != nothing
            top = popfirst!(arr)
            pushfirst!(arr, "")
            pushfirst!(arr, top)
        end
        tmp = table( [molname], [arr[2]], [category], [parse(Int64, arr[1])], [join(arr[3:end], " ")], names=cnames )
        #append!( rows(coeffs), rows(tmp) )
        coeffs = merge(coeffs, tmp)
    end
    
    Template(atoms, topologies, coeffs)
end

function addend(template, offset, nmol)
    new = Template()

    natom = select(template.atoms, :atomid) |> maximum
    adict = Dict()
    #:coeff注のatom typeの変換
    for imol in 1:nmol
        # !!! copy(template.atoms)とするとxのポインタ列をコピーするのでNG !!!
        add = copy(template).atoms
        select(add, :molid) .+= offset.molid + (imol-1)
        for i in 1:length(add)
            row = add[i]
            oldid = row[:atomid]
            column(add, :atomid)[i] += offset.atomid + (imol-1)*natom
            push!(adict, (row[:molid], oldid) => column(add, :atomid)[i])

            oldid = row[:atomtype]
            column(add, :atomtype)[i] += offset.atomtype
        end
        new.atoms = merge(new.atoms, add)
    end

    cdict = Dict()
    add = copy(template).coeffs
    for cat in select(template.coeffs, :category) |> unique, i in 1:length(add)
        row = add[i]
        if row[:category]==cat
            oldid = row[:cid]
            column(add, :cid)[i] += offset.cid[cat]
            push!(cdict, (cat, oldid)=>column(add, :cid)[i] )
        end
    end
    new.coeffs = merge(new.coeffs, add)

    categories = select(template.topologies, :category) |> unique
    ntop = Dict( cat => filter(t->t[:category]==cat, template.topologies) |> t->maximum(select(t, :tid)) for cat in categories )
    for imol in 1:nmol
        add = copy(template).topologies
        select(add, :molid)  .+= offset.molid + (imol-1)
        for cat in select(add, :category) |> unique,  i in 1:length(add)
            row = add[i]
            if row[:category] == cat
                cid_cat = cat2cat[cat]
                column(add, :cid)[i]    = cdict[ (cid_cat, row[:cid]) ]
                column(add, :atomid)[i] = [ adict[ (row[:molid], row[:atomid][j]) ] for j in 1:length(row[:atomid]) ]
                column(add, :tid)[i]   += offset.tid[cat] + (imol-1)*ntop[cat]
            end
        end
        new.topologies = merge(new.topologies, add)
    end

    return new
end

function myfilter(t, conditions...)
    # filterをかけるたびにコピーを生成
    tt = copy(t)
    for c in conditions
        key, f = c[1], c[2]
        tt = filter(f, tt, select=key)
    end
    tt
end

function alter!(t, mapping, conditions...)
    for row in t
        alltrue = true
        for c in conditions
            key, f = c[1], c[2]
            alltrue = alltrue && f(row[key])
        end
        if alltrue
            key = mapping[1]
            row[key] = mapping[2](row[key])
        end
    end
end








end #module