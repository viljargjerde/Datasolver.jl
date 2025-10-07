
using Datasolver, Revise, LinearAlgebra, Test, Plots

include("examples/simple_examples/addFuncs.jl")

#---------------------------------------------------------------------
### benchmark the Kanno truss with linear solution computed with 
### https://valdivia.staff.jade-hs.de/fachwerk_en.html
### solving with mixed 3-field formulation ######


# material
A = 2000/1e6        # [m²]
bar_E = 1.622e+07   # [Pa]

α = 1.0         # 0: linear strain     1: nonlinear strain

# load steps and factors
Fnodal = -40.0       # [N]
num_load_steps = 10
loadFac = LinRange(0.0,1.0,num_load_steps+1)

# geometry
node_vector = [
    [0,     0],
    [3.6,   0],
    [2*3.6, 0],
    [0,     3.6],
    [3.6,   3.6],
    [2*3.6, 3.6]
]  

constrained_dofs = [
    (1, 1),
    (1, 2),
    (4, 1),
    (4, 2)
]

connections = [
    (1, 2),
    (2, 3),
    (1, 5),
    (2, 4),
    (2, 5),
    (2, 6),
    (3, 5),
    (3, 6),
    (4, 5),
    (5, 6)
]


# force
force = zeros(2 * length(node_vector))
force[4] = Fnodal   # [N]   - downward force at node 2
force[6] = Fnodal   # [N]   - downward force at node 3

# truss problem
initProblem = TrussProblem(
    A,
    force,
    connections,
    α,
    constrained_dofs,
    node_vector = node_vector,
    num_quad_pts = 2,
)


results = solveANLP(
    initProblem=initProblem,
    constrained_dofs_global=constrained_dofs,
    externalForce=force,
    num_load_steps=num_load_steps,
    loadFac=Vector(loadFac),
    YoungModulus = bar_E,
    NR_max_iter=20,
    qrFactorized = false
);



# extract solution
uh = results.u[end]
uxh = uh[1:2:end]
uyh = uh[2:2:end]

eh = results.e[end]
sh = results.s[end]


# linear solution from online truss caculator
uRef = [
 	 0.000 	    0.000;
 	-9.084e-4 	-2.222e-3;
 	-1.174e-3 	-4.858e-3;
 	 0.000 	    0.000;
 	 8.672e-4 	-2.065e-3;
 	 1.045e-3 	-4.680e-3]

sRef = [-4093
		-1198
		-2697
		2960
		709.8
		-1135
		1694
		802.5
		3907
		802.5]

eRef = sRef ./ bar_E


# plot deformed structure
sc = 5e1
plot(0,0,dpi=150,size=(800,600),framestyle=:box)

for i in 1:length(connections)
    i1,i2 = connections[i]
    xN = [node_vector[i1][1],node_vector[i2][1]]
    yN = [node_vector[i1][2],node_vector[i2][2]]

    plot!(xN,yN,linewidth=2,linecolor=:black)

    xN = [node_vector[i1][1]+uxh[i1]*sc,node_vector[i2][1]+uxh[i2]*sc]
    yN = [node_vector[i1][2]+uyh[i1]*sc,node_vector[i2][2]+uyh[i2]*sc]

    plot!(xN,yN,linewidth=2,linecolor=:crimson)

    xN = [node_vector[i1][1]+uRef[i1,1]*sc,node_vector[i2][1]+uRef[i2,1]*sc]
    yN = [node_vector[i1][2]+uRef[i1,2]*sc,node_vector[i2][2]+uRef[i2,2]*sc]

    plot!(xN,yN,linewidth=2,linecolor=:royalblue)
end
plot!(legend=false,dpi=150, framestyle=:box, size=(800,600), xlabel="x", ylabel="y", tickfont=font(16), guidefont=font(16),legendfont=font(18))




# plot stress
plot(1:11, [sh[1];sh], linewidth=2,linetype=:steppre,label="ddcm",linecolor=:crimson)
plot!(1:11, [sRef[1];sRef], linewidth=2,linetype=:steppre, label="lin. ref. sol.",linecolor=:royalblue)

plot!(dpi=150, framestyle=:box, size=(800,600), xticks=(1.5:1:11,["1","2","3","4","5","6","7","8","9","10"]), xlabel="Element number", ylabel="Axial stress", tickfont=font(16), guidefont=font(16),legendfont=font(18), legend=:bottomright)



#---------------------------------------------------------------------





#---------------------------------------------------------------------
### load deflection curve

num_load_steps = 150
loadFac = LinRange(0.0,15.0,num_load_steps+1)


results = solveANLP(
    initProblem=initProblem,
    constrained_dofs_global=constrained_dofs,
    externalForce=force,
    num_load_steps=num_load_steps,
    loadFac=Vector(loadFac),
    YoungModulus = bar_E,
    NR_max_iter = 100,
    qrFactorized = false
);

# extract displacement of nodes under force and stress in elements 1+9
u2 = zeros(num_load_steps,2)
u3 = zeros(num_load_steps,2)

s1, s9 = zeros(num_load_steps), zeros(num_load_steps)

for i = 1:num_load_steps
    uh = results.u[i]

    u2[i,:] = uh[3:4]
    u3[i,:] = uh[5:6]
    s1[i] = results.s[i][1]
    s9[i] = results.s[i][9]
end


ll = loadFac[2:end]
plot(ll,u2[:,2],linewidth=2,linecolor=:crimson,label="uy, node 2")
plot!(ll,u3[:,2],linewidth=2,linecolor=:magenta,label="uy, node 3")

plot!(dpi=150, framestyle=:box, size=(800,600), xlabel="load factor", ylabel="uh", tickfont=font(16), guidefont=font(16),legendfont=font(18), legend=:bottomleft)

plot!(ll, -2.222e-3 .* ll,linewidth=1,linecolor=:black, linestyle=:dash, label="linear curve")
plot!(ll, -4.858e-3 .* ll,linewidth=1,linecolor=:black, linestyle=:dash, label=false)




plot(ll,s1,linewidth=2,linecolor=:crimson,label="ele 1")
plot!(ll,s9,linewidth=2,linecolor=:magenta,label="ele 9")

plot!(dpi=150, framestyle=:box, size=(800,600), xlabel="load factor", ylabel="Axial stress", tickfont=font(16), guidefont=font(16),legendfont=font(18))


plot!(ll, 4093 .* ll,linewidth=1,linecolor=:black, linestyle=:dash, label="lin. s1")
plot!(ll, -3907 .* ll,linewidth=1,linecolor=:black, linestyle=:dash, label="lin. s9")



#---------------------------------------------------------------------