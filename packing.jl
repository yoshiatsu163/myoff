module PackingMolecules
#export init_position

using JuliaDB, Statistics, Random

# :molname, :molid, :atomid, :atomtype, :q, :x
function init_position(atoms)
    # bounding box作成と同時に分子の重心を原点に移動
    bounds = [get_bounds(atoms, molname)[1] for molname in select(atoms, :molname) |> unique]
    bounds = [bounds[i][j] for i in 1:length(bounds), j in 1:3]
    bounds = [maximum(bounds[:,i]) for i in 1:3] .+ [1.0, 1.0, 1.0]

    idx = [findfirst(==(name), select(atoms, :molname)) for name in select(atoms, :molname) |> unique]
    push!(idx, length(atoms))
    molids = select(atoms, :molid)
    nmol = [molids[idx[i]]-molids[idx[i-1]] for i in 2:length(idx)]
    nmol[end] += 1

    @assert maximum(molids)==sum(nmol)
    npart = sum(nmol)^(1/3) |> ceil |> Int
    box = [bounds[i]*npart for i in 1:3]

    centers = [get_bounds(atoms, molname)[2] for molname in select(atoms, :molname) |> unique]
    cnt = 1
    for name in select(atoms, :molname) |> unique
        for i in 1:length(atoms)
            if atoms[i][:molname] == name
                column(atoms, :x)[i] .-= centers[cnt]
            end
        end
        cnt+=1
    end

    #move = [randperm(npart) randperm(npart) randperm(npart)]
    #cnt = 1
    #for mid in select(atoms, :molid) |> unique
    #    transport = [ (move[cnt,1]-0.5)*bounds[1], (move[cnt,2]-0.5)*bounds[2], (move[cnt,3]-0.5)*bounds[3] ]
    #    for iatom in 1:length(atoms)
    #        if atoms[iatom][:molid] == mid
    #            column(atoms, :x)[iatom] .+= transport
    #        end
    #    end
    #    cnt+=1
    #end

    #
    isvoid = vcat([false for i in 1:sum(nmol)],
                    [true  for i in sum(nmol)+1:npart^3]) |> shuffle!
    delta = [randperm(npart), randperm(npart), randperm(npart)]
    imol = 1
    for i in 1:npart, j in 1:npart, k in 1:npart
        if isvoid[ (i-1)*npart^2 + (j-1)*npart + k ]
            continue
        end
        for iatom in 1:length(atoms)
            if atoms[iatom][:molid] == imol
                column(atoms, :x)[iatom] .+= [ (i-0.5)*bounds[1], (j-0.5)*bounds[2], (k-0.5)*bounds[3] ]
            end
        end
        imol += 1
    end
    @assert imol == sum(nmol)+1

    return atoms, box
end

function get_bounds(atoms, molname)
    mols = filter(t->t[:molname]==molname, atoms)
    mol  = filter(t->t[:molid]==select(mols,:molid)[1], mols)
    p = column(mol, :x)

    p   = [p[i][j] for i in 1:length(p), j in 1:3]
    center = [0.5*(maximum(p[:,i])+minimum(p[:,i])) for i in 1:3]
    bounds = [maximum(p[:,i])-minimum(p[:,i]) for i in 1:3]

    (bounds, center)
end


    
end #module



