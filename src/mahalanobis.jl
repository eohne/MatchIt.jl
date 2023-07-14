

using Distances, CSV, RCall

data = CSV.File("Example Data/lalonde.csv") |> DataFrame;
treat = data[isequal.(data.treat, 1),[:age, :educ, :married,:nodegree,:re74]] |> Matrix
control = data[isequal.(data.treat, 0),[:age, :educ, :married,:nodegree,:re74]] |> Matrix



sqmahalanobis.((treat[1,:],),control,(cov(control),))

distance = Float64[]
index = Int[]
for j in 1:size(treat,1)
    dist = Float64[]
    for i in 1:size(control,1)
        push!(dist, sqmahalanobis(treat[j,:],control[i,:],cov(control)))
    end 
    temp_dist, idx = findmin(dist)
    push!(distance, temp_dist)
    push!(index,idx)
end
index
controls = control[index,:]

for i in 1:size(controls,2)
    print("\n$(mean(control[i]))\t\t $(mean(controls[i]))\t\t $(mean(treat[i]))")
end



# To Figure out: 
# 1. which covariance to take!
# 2. what is X and what is Y!
# 3. make more efficient (double loop is not efficient)
# 4. Compare to MatchIt




# Figure out:
# without replacement search:
#### 




