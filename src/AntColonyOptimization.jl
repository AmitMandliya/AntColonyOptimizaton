module AntColonyOptimization

using Random
using Printf
using Dates

# convert the msa to tsp graph. Takes sequences and input and returns the distance matrix
function MSA_to_TSP(sequences)
  node = []
  for i in 1:length(sequences)
    push!(node,sequences[i])
  end

  graph = Array{Float64, 2}(undef, length(node), length(node))
  max_score = 0
  for i in 1:length(node)
    graph[i,i] = 0
    for j in i+1:length(node)
      score, align1, align2 = produce_global_alignment(node[i],node[j]) 
      graph[i,j] = score
      graph[j,i] = score
      if max_score < score 
        max_score = score
      end
    end
  end
  # convert score to distance, more the score lesser the distance
  for i in 1:length(node)
    for j in i+1:length(node)
      graph[i,j] = max_score - graph[i,j] + 1
      graph[j,i] = max_score - graph[j,i] + 1
    end
  end    
  return graph
end

# perform the needleman-wunsch on two sequences v and w, returns the score and algined sequences
function produce_global_alignment(v, w, match_penalty=1, mismatch_penalty=-1, deletion_penalty=-1)
    rows = length(v)+1
    cols = length(w)+1
    dp = zeros(Float64, rows, cols)
    direction = zeros(Float64, rows, cols)
    for i in 1:(rows)
        dp[i,1] = (i-1) * deletion_penalty
        direction[i,1] = 2
    end
    for i in 1:(cols)
        dp[1,i] = (i-1) * deletion_penalty
        direction[1,i] = 3
    end
    for i in 2:(rows)
        for j in 2:(cols)
            if v[i-1] == w[j-1]
                ms = dp[i-1,j-1] + match_penalty
            else
                ms = dp[i-1,j-1] + mismatch_penalty
            end
            score = [ms, dp[i-1,j] + deletion_penalty, dp[i,j-1] + deletion_penalty]
            max_element = findmax(score)
            dp[i,j] = max_element[1]
            direction[i,j] = max_element[2]
        end
    end
    i = rows
    j = cols
    aligned_v = []
    aligned_w = []
    while(i > 1 || j > 1)
        dir = direction[i,j]
        if dir == 1
            i = i-1
            j = j-1
            push!(aligned_v, v[i])
            push!(aligned_w, w[j])
        end
        if dir == 2
            i=i-1
            push!(aligned_v, v[i])
            push!(aligned_w, "-")
        end
        if dir == 3
            j = j-1
            push!(aligned_v, "-")
            push!(aligned_w, w[j])
        end
    end
    return (dp[rows,cols], join(reverse(aligned_v)), join(reverse(aligned_w)))
end

# returns the ant colony with cities initialized to each ant
function create_colony(num_ants, num_nodes)
    colony = []
    for i in 1 : num_ants
       push!(colony, Dict("path"=>[rand(1:num_nodes)], "distance" => 0))
    end
    return colony
end

# create the pheror matrix
function create_pheror_matrix(num_nodes)
    pheromone = zeros(Float64, num_nodes, num_nodes)
    for i in 1: num_nodes
        for j in 1: num_nodes
            pheromone[i,j] = 1/num_nodes
        end
    end
    return pheromone
end

# intialize the probability matrix
function calculate_proba(num_nodes, pheromone, distance_matrix, alpha, beta)
    probability = zeros(Float64, num_nodes, num_nodes)
    for i in 1: num_nodes
        for j in 1: num_nodes
            probability[i,j] = (pheromone[i,j]^alpha) * (distance_matrix[i,j]^-beta)
            probability[j,i] = probability[i,j]
        end
    end
    return probability
end

# calculate probability of each ant going to all other nodes from current node
function calculate_proba_ant(pheromone, distance_matrix, unvisited_nodes, current_node, proba, alpha, beta)
  sigma = 0.0
  for unvisited_node in unvisited_nodes
    sigma += (pheromone[current_node,unvisited_node]^alpha) * (distance_matrix[current_node,unvisited_node]^-beta)
  end
  proba_ant = proba[current_node,:]/sigma
  return proba_ant
end

# finds the best path out of all ants
function find_best_path(n_ants, colony)
  bpath = []
  best_distance = Inf32
  typeof(best_distance)
  idx_best = 0
  for i=1: n_ants
    if colony[i]["distance"] < best_distance
      best_distance = colony[i]["distance"]
      bpath = colony[i]["path"]
      idx_best = i
    end
  end
  best_path = Dict("path"=> bpath, "distance"=> best_distance, "ant"=> idx_best)
  return best_path
end

# updates the pheror matrix
function update_pheror_matrix(num_nodes, n_ants, pheromone, distance_matrix, colony, Q, decay)
  depositpher = 0.0
  for i=1: n_ants
    ant = i
    for j= 1:(length(colony[ant]["path"])-1)
      src = colony[ant]["path"][j]
      dest = colony[ant]["path"][j+1]
      pheromone[src,dest] += Q/colony[i]["distance"]
    end
    depositpher += Q/colony[i]["distance"]
    for i= 1:num_nodes
      for j= 1:num_nodes
        pheromone[i,j] = (1-decay)*pheromone[i,j]*depositpher
        pheromone[j,i] = pheromone[i,j]
      end
    end
  end
  return pheromone
end

# calculate distance for one ant
function calculateDist_ant(ant, colony, distmatrix)
  dist = 0
  path = colony[ant]["path"]
  for i= 1:length(path)-1
    dist += distmatrix[path[i],path[i+1]]
  end
  return dist
end

# run each ant from source to destination
function traverse(ant, num_nodes, colony, pheromone, distance_matrix, proba, alpha, beta)
    unvisited = collect(1:num_nodes)
    current = colony[ant]["path"][1]
    deleteat!(unvisited, findfirst(isequal(current), unvisited))
    for j in 1: num_nodes-1
        if length(unvisited) > 1
            ant_probability = calculate_proba_ant(pheromone, distance_matrix, unvisited, current, proba, alpha, beta)
            prob = map((x) -> ant_probability[x] , unvisited)
            current = unvisited[findmax(prob)[2]]
            deleteat!(unvisited, findfirst(isequal(current), unvisited))
            push!(colony[ant]["path"], current)       
        else
            push!(colony[ant]["path"], unvisited[1])
        end
    end
    colony[ant]["distance"] = calculateDist_ant(ant, colony, distance_matrix)
end

# driver function for aco
function aco_main(num_ants, num_nodes, distance_matrix, iterations, Q, decay, alpha, beta)
    pheromone = create_pheror_matrix(num_nodes)
    gbpath = Dict()
    for i= 1: iterations
        colony = create_colony(num_ants, num_nodes)
        probability = calculate_proba(num_nodes, pheromone, distance_matrix, alpha, beta)
        for ant in 1: num_ants
            traverse(ant, num_nodes, colony, pheromone, distance_matrix, probability, alpha, beta)
        end
        # complete update_pheromone_matrix
        pheromone = update_pheror_matrix(num_nodes, num_ants, pheromone, distance_matrix, colony, Q, decay)
        #complete find best path fucntion
        best_path = find_best_path(num_ants, colony)
        
        bpath = best_path
        if i == 1
            gbpath = bpath
        else
            if bpath["distance"] < gbpath["distance"]
                gbpath = bpath
            end
        end
    end 
    # return the best path
    return gbpath["path"]
end

function find_gap_indices(A, alignedA)
    i = 1
    j = 1
    pointer = []
    while j <= length(alignedA)
        if alignedA[j] == '-' && (i > length(A) || A[i] != '-')
            push!(pointer, j)
            j += 1
        else
            j += 1
            i += 1
        end
    end
    return pointer
end

function insert_gaps(S,gap_indices_for_A)
    copy_of_S = S
    if length(gap_indices_for_A) > 0 && length(gap_indices_for_A) > 0
        gap_indices_for_A = sort(gap_indices_for_A)
        for i in gap_indices_for_A
            copy_of_S = string(string(copy_of_S[1:i-1],'-'),copy_of_S[i:end])
        end
    end
    return copy_of_S
end

# this function takes two params: sequences is the original array of input sequences, order is the array output of tsp algorithm of the order of sequences
# notice that the index in both the params are relative to each other
function align_output_sequences(sequences, order) 
    ordered_sequences = Array{String,1}(undef,0)
    for i=1:length(order)
        push!(ordered_sequences,sequences[order[i]])
    end
    aligned_sequences = Array{String,1}(undef,0)
    for i=1:length(ordered_sequences)
        push!(aligned_sequences,ordered_sequences[i])
    end
    for i=1:length(aligned_sequences)-1
        A = aligned_sequences[i]
        B = aligned_sequences[i+1]
        score, alignedA, alignedB = produce_global_alignment(A,B)
        gap_indices_for_A = find_gap_indices(A, alignedA)
        # go to all predecessors of A and insert the gaps at same place
        for j = 1:i-1
            S = aligned_sequences[j]
            newly_alinged_S = insert_gaps(S,gap_indices_for_A)
            aligned_sequences[j] = newly_alinged_S
        end
        aligned_sequences[i] = alignedA
        aligned_sequences[i+1] = alignedB
    end  
    ordered_aligned_sequences = Array{String,1}(undef,length(order))
    for i=1:length(order)
        ordered_aligned_sequences[order[i]] = aligned_sequences[i]
    end
    return ordered_aligned_sequences
end

# calculate the final score
function calc_sum_pair(sequences) 
    t = length(sequences)
    k = length(sequences[1])
    score = 0
    for i=1:t
        A = sequences[i]
        for j=i+1:t
            B = sequences[j]
            for idx = 1:k
                if A[idx] == B[idx] && A[idx] != '-'
                    score += 1
                end
            end
        end
    end
    return score
end

# calculate the final score to include gap penalty
function calc_sum_pair_include_gap_penalty(sequences) 
    t = length(sequences)
    k = length(sequences[1])
    score = 0
    for i=1:t
        A = sequences[i]
        for j=i+1:t
            B = sequences[j]
            for idx = 1:k
                # Add 1 for match
                if A[idx] == B[idx] && A[idx] != '-'
                    score += 1
                # subtract 1 for mismatch
                elseif A[idx] != B[idx] && A[idx] != '-' && B[idx] != '-'
                    score -= 1
                # subtract 1 for gap
                elseif A[idx] == '-' ||  B[idx] == '-' 
                    score -= 1
                end
            end
        end
    end
    return score
end

# get sequences from the input file
function get_sequences_from_file(file_name)
    sequences = []
    flag_first_arrow = true
    open(file_name) do f
        sequence = ""
        for line in eachline(f)
            if startswith(line, '>') 
                if flag_first_arrow
                    flag_first_arrow = false
                    continue
                end
                push!(sequences, sequence)
                sequence = ""
            else
                sequence *= line
            end
        end
        if length(sequence) > 0
            push!(sequences, sequence)
        end
    end
    return sequences
end

# write output to the file
function write_to_file(output_filename, aligned_sequences)
    touch(output_filename)
    f = open(output_filename, "w") 
    for i=1:length(aligned_sequences)
        seq_num = string(">s",i)
        write(f,string(seq_num,"\n"))
        write(f,aligned_sequences[i])
        write(f,"\n")
    end
    close(f) 
end

# this is the driver function which take the fasta file as input and writes the output to the output.txt file.
# the default input to this function is input.fasta which should be present in the same folder as this code.
function driver(filename="input1.txt",output_filename = "output1.txt")
    # get the sequences from input file
    input_sequences = get_sequences_from_file(filename)
    # convert the get input_sequences to TSP problem and get the distance_matrix for TSP
    input_distance_matrix = MSA_to_TSP(input_sequences)
    # configuration for ACO
    num_sequences = length(input_sequences)
    num_ants = maximum([100,5*num_sequences])
    iterations = 30
    Q = 0.6
    decay = 0.6
    alpha = 1
    beta = 1
    # run the ant colony optimization, get the order of node traversal
    order_of_alignment = aco_main(num_ants,num_sequences,input_distance_matrix,iterations,Q,decay,alpha,beta)
    # align output sequences in the order identifiedgi by above algorithm
    aligned_sequences = align_output_sequences(input_sequences,order_of_alignment)
    # get the score of the aligned_sequences
    final_score = calc_sum_pair(aligned_sequences)
    # save the output in fasta format in output1.txt
    write_to_file(output_filename, aligned_sequences)
    # print the aligned sequences
    println("Score of aligned sequences = ",final_score)
    println("The aligned sequences are saved in output file. Printing them here as well:")
    for i=1:length(aligned_sequences)
        println(aligned_sequences[i])
    end  
    return final_score, aligned_sequences
end

end