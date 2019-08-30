library(ggplot2)

cpu_dense = read.csv("../timings/cpu_dense_times.csv")
cpu_sparse = read.csv("../timings/cpu_sparse_times.csv")
gpu_Tesla_dense = read.csv("../timings/gpu_Tesla_dense_times.csv")
gpu_Tesla_sparse = read.csv("../timings/gpu_Tesla_sparse_times.csv")
gpu_Tesla_dense = read.csv("../timings/gpu_Tesla_dnsspr_times.csv")
gpu_Tesla_femses = read.csv("../timings/gpu_Tesla_femses_times.csv")

#cpu_sparse[5][cpu_sparse[1] == 4] / gpu_Tesla_sparse[5][gpu_Tesla_sparse[1] == 4]

get_cpu_speedups <- function(gpu, cpu){
  prob_sizes = as.numeric(unlist(unique(cpu[1])))
  speedups = gpu
  for(i in prob_sizes){
    for(j in 5:12){
      speedups[j][speedups[1] == i] = cpu[j][cpu[1] == i] / gpu[j][gpu[1] == i]
    }
  }
  
  return(speedups)
}

get_gpu_speedups <- function(gpu1, gpu2){
  prob_sizes = as.numeric(unlist(unique(gpu2[1])))
  speedups = gpu2
  for(i in prob_sizes){
    for(j in 5:12){
      speedups[j][speedups[1] == i] = gpu2[j][gpu2[1] == i] / gpu1[j][gpu1[1] == i]
    }
  }
  
  return(speedups)
}

###############
#get_gpu_speedups <- function(new)

# ggplot(cpu_dense, aes(n, total)) + geom_point()
# ggplot(cpu_sparse, aes(n, total)) + geom_smooth() + geom_point()

#ggplot(gpu_Tesla_sparse, aes(block_size_X, total, colour = as.factor(n+1))) + geom_point() + geom_smooth()
# 
# 
# prob_sizes = as.numeric(unlist(unique(gpu_Tesla_dense[1])))
# prob_sizes
# for(i in prob_sizes){
#   print(i)
#   print(max(gpu_Tesla_dense[5][gpu_Tesla_dense[1] == i] / gpu_Tesla_sparse[5][gpu_Tesla_sparse[1] == i]))
# }
###############

speedups = get_cpu_speedups(gpu_Tesla_sparse, cpu_sparse)
ggplot(speedups, aes(block_size_X, log10(total), colour = as.factor(n+1))) + geom_point() + geom_smooth()
ggplot(speedups, aes(block_size_X, log10(assembly), colour = as.factor(n+1))) + geom_point() + geom_smooth()
ggplot(speedups, aes(block_size_X, log10(elem_mats), colour = as.factor(n+1))) + geom_point() + geom_smooth()

speedups = get_cpu_speedups(gpu_Tesla_dense, cpu_dense)
ggplot(speedups, aes(block_size_X, log10(total), colour = as.factor(n+1))) + geom_point() + geom_smooth()
ggplot(speedups, aes(block_size_X, log10(assembly), colour = as.factor(n+1))) + geom_point() + geom_smooth()
ggplot(speedups, aes(block_size_X, log10(elem_mats), colour = as.factor(n+1))) + geom_point() + geom_smooth()

speedups = get_gpu_speedups(gpu_Tesla_sparse, gpu_Tesla_dense)
ggplot(speedups, aes(block_size_X, log10(total), colour = as.factor(n+1))) + geom_point() + geom_smooth()
### use these two plots to demonstrate eratticity of read/write of timings array
ggplot(speedups, aes(block_size_X, log10(assembly), colour = as.factor(n+1))) + geom_point() + geom_smooth()
ggplot(speedups, aes(block_size_X, log10(elem_mats), colour = as.factor(n+1))) + geom_point() + geom_smooth()


ggplot(speedups, aes(block_size_X, log10(assembly), colour = as.factor(reconfig))) + geom_point() + geom_smooth()

##############

gpu_Tesla_sparse = read.csv("~/Desktop/timings2/gpu_Tesla_sparse_times.csv")

ggplot(subset(gpu_Tesla_sparse,reconfig==0), aes(block_size_X, total, colour = as.factor(reconfig))) + geom_point() + geom_smooth()
ggplot(gpu_Tesla_sparse, aes(block_size_X, assembly, colour = as.factor(reconfig))) + geom_point() + geom_smooth()
ggplot(gpu_Tesla_sparse, aes(block_size_X, elem_mats, colour = as.factor(reconfig))) + geom_point() + geom_smooth()
ggplot(gpu_Tesla_sparse, aes(block_size_X, elems_p_assemb, colour = as.factor(reconfig))) + geom_point() + geom_smooth()
ggplot(gpu_Tesla_sparse, aes(block_size_X, sparsity.scan)) + geom_point() + geom_smooth()
