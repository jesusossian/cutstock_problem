using JuMP
using GLPK
using Gurobi
using SparseArrays
using LinearAlgebra

function solve_pricing(dual_rest, maxwidth, widths, demand)
    reduced_costs = dual_rest
    n = length(reduced_costs)
    
    # The actual pricing model.
    submodel = Model(Gurobi.Optimizer)
    set_silent(submodel)
    @variable(submodel, xs[1:n] >= 0, Int)
    @constraint(submodel, sum(xs .* widths) <= maxwidth)
    @objective(submodel, Min, 1 - sum(xs .* reduced_costs))
    
    println(submodel)
        
    # Solve the sub problem
    optimize!(submodel)
    
    new_pattern = round.(Int, value.(xs))
    
    net_cost = objective_value(submodel)
    println("Reduced cost: ", net_cost)
    println("Reduced cost column= ", new_pattern)
       
    # If the net cost of this new pattern is nonnegative, no more patterns to add.
    return net_cost >= 0 ? nothing : new_pattern

end # end function solve_pricing

function example_cutting_stock(; max_gen_cols::Int = 5_000)

    maxwidth = 100
    m = 5    
    A = Array{Int64}(I,m,m)
       
    println("maxwith = $(maxwidth), m = $(m)")
    
    widths = [22; 42; 52; 53; 78]
    demand = [45; 38; 25; 11; 12]
    
    n = length(widths)
    ncols = length(widths) # initial number of variables
    
    println("n = $(n), ncols = $(ncols)")
    
    readline()
           
    # Write the master problem with this "reduced" set of patterns.
    # Not yet integer variables: otherwise, the dual values may make no sense
    # (actually, Guobi(or GLPK) will yell at you if you're trying to get duals for
    # integer problems).

    println("# First master problem")
    master = Model(Gurobi.Optimizer) # modelo for the master problem
    set_silent(master)
    @variable(master, x[1:ncols] >= 0) # setting variables
    @objective(master, Min, sum(1 * x[j] for j=1:ncols)) # setting objective function
    @constraint(master, rest[i=1:m], sum(A[i,j] * x[j] for  j=1:ncols) == demand[i])
    println(master)
    
    # First solve of the master problem.
    optimize!(master)

    println("Solution of the master problem: ", value.(x))
    println("Objective value of the master problem: ", objective_value(master))
    
    if termination_status(master) != MOI.OPTIMAL
        warn("Master not optimal ($ncols patterns so far)")
    end
    
    readline()
    
    # Then, generate new patterns, based on the dual information.
    while ncols - n <=  max_gen_cols ## Generate at most max_gen_cols columns.
                              
        if ! has_duals(master)
            break
        end
        
        println("\n# The actual pricing problem ")
        new_pattern = solve_pricing(dual.(rest), maxwidth, widths, demand)
        
        readline()
        
        # No new pattern to add to the formulation: done!
        if new_pattern === nothing
            break
        end
        
        # Otherwise, add the new pattern to the master problem, recompute the
        # duals, and go on waltzing one more time with the pricing problem.
        ncols += 1
        A = hcat(A, new_pattern)
        
        #println(A)
                
        # One new variable.
        new_var = @variable(master, [ncols], base_name = "x", lower_bound = 0)
        push!(x, new_var[ncols])
        
        # Update the objective function.
        set_objective_coefficient(master, x[ncols], 1.0)
        
        # Update the constraint number j if the new pattern impacts this production.
        for j = 1:n
            if new_pattern[j] > 0
                set_normalized_coefficient(rest[j], new_var[ncols], new_pattern[j])
            end
        end
        
 
        println("# The actual master problem: pattern $(ncols)\n")
        println(master)
                
        # Solve the new master problem to update the dual variables.
        optimize!(master)
 
        println("Solution of the master problem = ", value.(x))
        println("Objective value of the master problem = ", objective_value(master))
    
        readline()
            
        if termination_status(master) != MOI.OPTIMAL
            @warn("Master not optimal ($ncols patterns so far)")
        end
        
    end
    
    # Finally, impose the master variables to be integer and resolve.
    # To be exact, at each node in the branch-and-bound tree, we would need to
    # restart the column generation process (just in case a new column would be
    # interesting to add). This way, we only get an upper bound (a feasible
    # solution).
    set_integer.(x)
    
    optimize!(master)
    
    if termination_status(master) != MOI.OPTIMAL
        @warn("Final master not optimal ($ncols patterns)")
        return
    end
    
    println("\n********** Final Solution **********\n")
    opt_master = objective_value(master)
    println("Objective value = $(opt_master)")
    println("Optimal solution = ", value.(x))
    
    println("\n# Number of each pattern: ")
    for i = 1:length(x)
        if value(x[i]) > 0.5
            println("$(round(Int, value(x[i]))) units of pattern $(i)")
        end
    end

    return

end

example_cutting_stock()
