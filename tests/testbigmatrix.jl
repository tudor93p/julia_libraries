sizes_bigH = [size(H(a...),1) for a in atoms_bigH] 

#function get_sector(bigmatrix,A,B)
#
#	i = argmax([A==a for a in atoms_bigH])
#
#	j = argmax([B==a for a in atoms_bigH])
#
#	i = i==1 ? 0 : sum(sizes_bigH[1:i-1])
#
#	j = j==1 ? 0 : sum(sizes_bigH[1:j-1])
#
#	n,m = size(H(A...,B...))
#
#	return bigmatrix[i+1:i+n,j+1:j+m]
#
#end


function couplingUG(atom)

	u = hcat([H(atom...,b...) for b in atoms_bigH]...)

	return (u,Graph.get_prop(g,node(atom...)...,:GF))

end

bruteforcegf = RGF(Energy,
									 hcat([vcat([H(a...,b...) for a in atoms_bigH]...) for b in atoms_bigH]...),couplingUG.(attached_atoms)...)



layeredgf = hcat([vcat([G(a...,b...) for a in atoms_bigH]...) for b in atoms_bigH]...)



println()
@show LA.norm(layeredgf - bruteforcegf) 

println()

