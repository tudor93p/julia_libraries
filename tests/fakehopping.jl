
#####################################333
d0=3
function fakehopping(ri,rj)

	d = LA.norm(ri.-rj)

	isapprox(d,0,atol=dist_tol) && return LA.diagm(0=>repeat([0.1+sum(ri)],d0))

	isapprox(d,1,atol=dist_tol) && return LA.diagm(0=>repeat([1.],d0)) +LA.diagm(1=>repeat((ri[1:1]-rj[1:1])*0.1im,2))


end



HoppMatr(A,B=A) = TBmodel.HoppingMatrix(A,B,Hopping=fakehopping,Nr_Orbitals=d0)


