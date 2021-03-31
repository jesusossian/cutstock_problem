using JuMP
using LinearAlgebra
using Gurobi

function solve_pricing(m, dual_rest, maxw, w)
  p = dual_rest
  n = length(p)
  
  # Describe the sub problem
  subprob = Model(Gurobi.Optimizer) # Model for the subproblem
  
  set_optimizer_attribute(subprob, "LogToConsole", 0)
  
  @variable(subprob, y[1:m] >= 0, Int )
  @objective(subprob, Min, 1 - sum(p[i] * y[i] for i=1:m) )
  @constraint(subprob, sum(w[i] * y[i] for i=1:m) <= maxw)
  print(subprob)
  
  # Solve the sub problem
  optimize!(subprob)

  red_cost = objective_value(subprob)

  if red_cost >= 0
    println("\nReduced cost = $(red_cost)")
    println("Solution of the master problem is optimal")
  else
    println("\nReduced cost = $(red_cost)")
    println("Cost reducing column = ", value.(y))
  end
 
  return red_cost, y

end


function cutting_stock()
  maxw = 100
  m = 5
  A = Array{Int64}(I,m,m)
  
  d = [45; 38; 25; 11; 12]
  w = [22; 42; 52; 53; 78]

  n = length(w)  
  ncols = length(w) # initial number of variables
  println("ncols = $(ncols)")
  
  t1 = time_ns()

  println("\nFirst master problem")
  master = Model(Gurobi.Optimizer) # model for the master problem
    
  set_optimizer_attribute(master, "LogToConsole", 0)
    
  set_silent(master)
  @variable(master, x[1:ncols] >= 0) # setting variables
  @objective(master, Min, sum(1 * x[j] for j=1:ncols)) # setting objective function
  @constraint(master, rest[i=1:m], sum(A[i,j] * x[j] for j=1:ncols) == d[i])

  println(master)
  
  # Solving the master problem with the initial data
  optimize!(master)
    
  opt_master = objective_value(master)
  println("\nSolution of the master problem: ", value.(x))
  println("Objective value of the master problem: ", opt_master)

  #Dual variables for the equality constraints
  p = zeros(m)
  for i=1:m
    p[i] = dual(rest[i]) #input data for the subproblem
  end 
  println("\nDual variables of master problem: ", p)

  readline()
        
  println("\nFirst subproblem")  
  red_cost, y = solve_pricing(m, p, maxw, w)
  
  k = 1
    
  while (red_cost < 0)

    # Newly added column to the A matrix
    Anew = Array{Float64}(undef, 0)
    for i=1:m
      push!(Anew, value.(y)[i])
    end
        
    readline()
     
    ncols += 1 
    
    println("\nNew master problem: pattern $(k) \n")
    new_var = @variable(master, [ncols], base_name = "x", lower_bound = 0)
    push!(x, new_var[ncols])
 
    set_objective_coefficient(master, x[ncols], 1.0)

    for i in 1:m
      if Anew[i] > 0
        set_normalized_coefficient(rest[i], new_var[ncols], Anew[i])
      end
    end
 
    println(master)
      
    optimize!(master)
    
    opt_master = objective_value(master)
    println("\nSolution of the master problem: ", value.(x))
    println("Objective value of the master problem: ", opt_master)

    #Dual variables for the equality constraints and store them in an array p
    for i in 1:m
      p[i] = dual(rest[i]) # These p[i] are the input data for the subproblem
    end 
    println("\nDual variables of master problem: ", p)

    readline()
 
    println("\nNew subproblem")      
    red_cost, y = solve_pricing(m, p, maxw, w)
      
    k = k+1

  end
    
  t2 = time_ns()
  elapsedtime = (t2 - t1)/1.0e9
  
  println("\nFinal solution:")
  println("Objective value end of master problem: ", objective_value(master))
  println("Current Solution end of master problem: ", value.(x))
  println("time = $(elapsedtime)")

end

cutting_stock()
