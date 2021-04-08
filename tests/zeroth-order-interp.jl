#!/home/tudor/apps/julia/julia-1.1.1/bin/julia
#

#import Algebra,Plots,Interpolations,Utils


x = sort(rand(10))


x2 = Utils.Rescale(sort(rand(300)),x)

y = rand(length(x))

#f(a) = Interpolations.interpolate(a, Interpolations.BSpline(Interpolations.Constant()))
Plots.scatter(x,y)


y2 = Algebra.Interp1D(x,y,0,x2)


Plots.plot!(x2,y2)
#
#
#
#y3 = f(y)(Utils.Rescale(x2,axes(x,1),x))
#
#Plots.plot!(x2,y3)




























