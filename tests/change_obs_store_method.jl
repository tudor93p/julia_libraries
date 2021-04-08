using JLD,DelimitedFiles, FileIO
import Utils

folder1="dat/"
folder2="jld/"

files = cd(readdir,folder1)

files1 = folder1.*files 

files2 = folder2.*files.*".jld"


for (f,f1,f2) in zip(files,files1,files2)

	if !occursin("Legend",f)
#
#		println(f)

#	 save(dest, "A", d["A"], "B", d["B"]) 
	
#	save(File(format"JLD",f2), Utils.Read_PhysObs(x->folder1*x,f))

	data = Utils.Read_PhysObs(x->folder1*x, f, "plain") 

#	@show keys(data)

	Utils.Write_NamesVals(x->folder2*x, "new", "jld")[1](f,data[f])


	data2 =  Utils.Read_PhysObs(x->folder2*x, f, "jld")


	for k in union(keys(data),keys(data2))

		d1,d2=data[k],data2[k]

		if d1 isa AbstractDict

			for q in union(keys(d1),keys(d2))

				println(maximum(abs.(d1[q]-d2[q])))

			end 

		else 

			println(maximum(abs.(d1-d2)))

		end 
	end 


#	@show round(filesize(f2)/filesize(f1)*100,digits=2)

#	println()
	println()

end 
end 
