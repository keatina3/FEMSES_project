library(ggplot2)
library(gridExtra)
library(scales)
library(reshape2)
library(wesanderson)

cpu_dense = read.csv("~/Documents/College/Project/src/timings/cpu_dense_times.csv")
cpu_sparse = read.csv("~/Documents/College/Project/src/timings/cpu_sparse_times.csv")
gpu_Tesla_dense = read.csv("~/Documents/College/Project/src/timings/gpu_Tesla_dense_times.csv")
gpu_Tesla_sparse = read.csv("~/Documents/College/Project/src/timings/gpu_Tesla_sparse_times.csv")
gpu_Tesla_dnsspr = read.csv("~/Documents/College/Project/src/timings/gpu_Tesla_dnsspr_times.csv")
gpu_Tesla_femses = read.csv("~/Documents/College/Project/src/timings/gpu_Tesla_femses_times.csv")

col_sprs = "#00AFBB"
col_dns = "#E7B800"
col_dnsspr = "#FC4E07"
col_femses = "#8A2BE2"

cpu_sparse[10] = cpu_sparse[8]  + cpu_sparse[9]
cpu_dense[10] = cpu_dense[8] + cpu_dense[9]

get_speedups <- function(gpu, cpu){
  prob_sizes = as.numeric(unlist(unique(cpu[1])))
  speedups = gpu
  for(i in prob_sizes){
    for(j in 5:12){
      speedups[j][speedups[1] == i] = cpu[j][cpu[1] == i] / gpu[j][gpu[1] == i]
    }
  }
  
  return(speedups)
}

get_averages <- function(df){
  avgs = data.frame(aggregate(df,by=list(df$n, df$reconfig!=0, df$block_size_X),data=df,FUN=mean))[4:18]
  
  return(avgs)
}

base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

cpu_sparse_avgs = get_averages(cpu_sparse)
cpu_dense_avgs = get_averages(cpu_dense)
gpu_Tesla_sparse_avgs = get_averages(gpu_Tesla_sparse)
gpu_Tesla_dense_avgs = get_averages(gpu_Tesla_dense)
gpu_Tesla_dnsspr_avgs = get_averages(gpu_Tesla_dnsspr) 
gpu_Tesla_femses_avgs = get_averages(gpu_Tesla_femses)


###### Total speedups over CPU vs problem_size ########

n_labeller <- function(variable,value){
  return(labels[value])
}

speedups = get_speedups(gpu_Tesla_sparse_avgs, cpu_sparse_avgs)
total_sparse_cpu_speedup_vs_n <- ggplot(subset(speedups,reconfig!=0), aes((n+1)*(n+1), total)) + 
                              geom_smooth(colour=col_sprs, formula = y~log(x)) + geom_point(colour=col_sprs)

speedups = get_speedups(gpu_Tesla_dense_avgs, cpu_dense_avgs)
total_dense_cpu_speedup_vs_n <- ggplot(subset(speedups,reconfig!=0), aes((n+1)*(n+1), total)) +
                            geom_smooth(colour = col_dns, formula = y~log(x)) + geom_point(colour = col_dns)

speedups = get_speedups(gpu_Tesla_dnsspr_avgs, cpu_sparse_avgs)
total_dnsspr_cpu_sparse_speedup_vs_n = ggplot(subset(speedups,reconfig!=0), aes((n+1)*(n+1), total))  +
                            geom_smooth(colour = col_dnsspr, formula = y~log(x)) + geom_point(colour=col_dnsspr)
speedups = get_speedups(gpu_Tesla_dnsspr_avgs, cpu_dense_avgs)
total_dnsspr_cpu_dense_speedup_vs_n <- ggplot(subset(speedups,reconfig!=0), aes((n+1)*(n+1), total)) +
                            geom_smooth(colour = col_dnsspr, formula = y~log(x)) + geom_point(colour=col_dnsspr)

data.m <- melt(subset(gpu_Tesla_spare_avgs), id.vars=c("n", "m", "total", "block_size_X", "reconfig", "sse", 
                                                          "iterations", "elems_p_assemb", "convert"))
ggplot(data.m, aes(fill=variable, y=value/total, x=as.factor((n+1)*(n+1)))) +
        geom_bar(position="fill", stat="identity") +
        labs(title="Proportion of Computation Time Taken By Each Step\nin Sparse Serial Process", 
        x ="Problem Size", y = NULL, fill = "Kernel") +
        scale_y_continuous(labels = scales::percent_format())

data.m <- melt(subset(gpu_Tesla_dense_avgs), id.vars=c("n", "m", "total", "block_size_X", "reconfig", "sse", 
                                                          "iterations", "elems_p_assemb", "convert"))
ggplot(data.m, aes(fill=variable, y=value/total, x=as.factor((n+1)*(n+1)))) +
        geom_bar(position="fill", stat="identity") +
        labs(title="Proportion of Computation Time Taken By Each Step\nin Sparse Serial Process", 
        x ="Problem Size", y = NULL, fill = "Kernel") +
        scale_y_continuous(labels = scales::percent_format())


## femses 
speedups = get_speedups(gpu_Tesla_femses_avgs, cpu_sparse_avgs)
total_femses_cpu_sparse_speedup_vs_n = ggplot(subset(speedups,reconfig!=0), aes((n+1)*(n+1), total)) +
                            geom_smooth(colour = col_sprs,formula = y~log(x)) + geom_point(colour=col_sprs)
    
speedups = get_speedups(gpu_Tesla_femses_avgs, cpu_dense_avgs)
total_femses_cpu_dense_speedup_vs_n <- ggplot(subset(speedups,reconfig!=0), aes((n+1)*(n+1), total)) + 
                            geom_smooth(colour = col_dns,formula = y~log(x)) + geom_point(colour=col_dns)

speedups = get_speedups(gpu_Tesla_femses_avgs, gpu_Tesla_sparse_avgs)
total_femses_gpu_sparse_speedup_vs_n = ggplot(subset(speedups,reconfig!=0), aes((n+1)*(n+1), total)) +
                            geom_smooth(colour = col_dnsspr, formula = y ~ log(x)) + geom_point(colour=col_dnsspr)
speedups = get_speedups(gpu_Tesla_femses_avgs, gpu_Tesla_dense_avgs)
total_femses_gpu_dense_speedup_vs_n <- ggplot(subset(speedups,reconfig!=0), aes((n+1)*(n+1), total)) +
                            geom_smooth(colour = col_femses, formula = y ~ log(x)) + geom_point(colour=col_femses)


## figures 

total_sparse_cpu_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                  labs(title="Speedup Over Serial Total Time vs Problem Size for Sparse Solver", x ="Degrees of Freedom", y = "Speedup")
total_dense_cpu_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                  labs(title="Speedup Over SerialTotal Time vs Problem Size for Dense Solver", x ="Degrees of Freedom", y = "Speedup")
total_dnsspr_cpu_sparse_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                  labs(title="Speedup Over Serial Sparse Solver Total Time vs Problem Size for\n Dense-Sparse Conversion Solver", 
                       x ="Degrees of Freedom", y = "Speedup")
total_dnsspr_cpu_dense_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                  labs(title="Speedup Over Serial Dense Solver Total Time vs Problem Size for\n Dense-Sparse Conversion Solver", 
                       x ="Degrees of Freedom", y = "Speedup")
total_femses_cpu_sparse_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                  labs(title="Speedup of FEMSES Total Time Over Serial Sparse Solver vs Problem Size", x ="Degrees of Freedom", y = "Speedup")
total_femses_cpu_dense_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                  labs(title="Speedup of FEMSES Total Time Over Serial Dense Solver vs Problem Size", x ="Degrees of Freedom", y = "Speedup")
total_femses_gpu_sparse_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                  labs(title="Speedup of FEMSES Total Time Over GPU Sparse Solver vs Problem Size", x ="Degrees of Freedom", y = "Speedup")
total_femses_gpu_dense_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                  labs(title="Speedup of FEMSES Total Time Over GPU Dense Solver vs Problem Size", x ="Degrees of Freedom", y = "Speedup")


###### Total speedup over CPU vs block_size #######

speedups = get_speedups(gpu_Tesla_sparse_avgs, cpu_sparse_avgs)
total_sparse_cpu_speedup_vs_b <- ggplot(subset(speedups,reconfig!=0), aes(block_size_X, total,
                                        colour=as.factor((n+1)*(n+1)))) + geom_smooth() + geom_point()
speedups = get_speedups(gpu_Tesla_dense_avgs, cpu_dense_avgs)
total_dense_cpu_speedup_vs_b <- ggplot(subset(speedups,reconfig!=0), aes(block_size_X, total,
                                        colour=as.factor((n+1)*(n+1)))) + geom_smooth() + geom_point()

## femses times
# separate times for M enabled/disabled??

speedups = get_speedups(gpu_Tesla_femses_avgs, cpu_sparse_avgs)
total_femses_cpu_sparse_speedup_vs_b <- ggplot(subset(speedups,reconfig!=0), aes(block_size_X, total,
                                        colour=as.factor((n+1)*(n+1)))) + geom_smooth(formula = y ~ log(x)) + geom_point()
speedups = get_speedups(gpu_Tesla_femses_avgs, cpu_sparse_avgs)
labels = unique(subset(speedups,n<200)[1])
labels = list((labels+1) **2 )
total_femses_cpu_sparse_speedup_vs_b_facet <- ggplot(subset(speedups, block_size_X > 10), 
                                        aes(block_size_X, total, colour = as.factor(reconfig))) +
                                        geom_smooth(colour=col_sprs) + geom_point(colour=col_sprs) +
                                        facet_wrap(~n, scales="free")

speedups = get_speedups(gpu_Tesla_femses_avgs, cpu_dense_avgs)
total_femses_cpu_dense_speedup_vs_b <- ggplot(subset(speedups,reconfig=!0), aes(block_size_X, total,
                                        colour=as.factor((n+1)*(n+1)))) + geom_smooth(formula = y ~ log(x)) + geom_point()
total_femses_cpu_dense_speedup_vs_b_facet <- ggplot(subset(speedups, reconfig!=0  & block_size_X > 10), aes(block_size_X, total)) +
                                        geom_smooth(colour=col_dns) + geom_point(colour=col_dns) +
                                        facet_wrap(~n, scales="free")
                        

## figures
total_sparse_cpu_speedup_vs_b + scale_y_log10(breaks = base_breaks()) +
                                      labs(title="Speedup Over Serial of Total Time vs Block Size for Sparse Solver", 
                                      x ="Block Size", y = "Speedup") +
                                      scale_color_discrete(name="DOF")
total_dense_cpu_speedup_vs_b  + scale_y_log10(breaks = base_breaks()) + 
                                    labs(title="Speedup Over Serial Total Time vs Block Size for Dense Solver", 
                                     x ="Block Size", y = "Speedup") +
                                    scale_color_discrete(name="DOF")
total_femses_cpu_sparse_speedup_vs_b + scale_y_log10(breaks = base_breaks()) + 
                                    labs(title="Speedup of FEMSES Total Time Over Serial Sparse Solver vs Block Size", 
                                    x ="Block Size", y = "Speedup") +
                                    scale_color_discrete(name="DOF")
total_femses_cpu_dense_speedup_vs_b + scale_y_log10(breaks = base_breaks()) +
                                    labs(title="Speedup of FEMSES Total Time Over Serial Dense Solver vs Block Size", 
                                    x ="Block Size", y = "Speedup") +
                                    scale_color_discrete(name="DOF")

total_femses_cpu_sparse_speedup_vs_b_facet + scale_y_log10(breaks = base_breaks()) + 
                              labs(title="Speedup of FEMSES Total Time Over Serial Sparse Solution\nvs Block Size for Individual Problem Sizes", 
                                x ="Block Size", y = "Speedup") +
                                scale_color_discrete(name="Mem Config")
total_femses_cpu_dense_speedup_vs_b_facet + scale_y_log10(breaks = base_breaks()) + 
                                labs(title="Speedup of FEMSES Total Time Over Serial Dense Solution\nvs Block Size for Individual Problem Sizes", 
                                x ="Block Size", y = "Speedup") +
                                scale_color_discrete(name="Mem Config")

###### Solver times vs problem size #####
# try get these on same plot???  for 4th plot#

speedups = get_speedups(gpu_Tesla_sparse_avgs, cpu_sparse_avgs)
solve_sparse_cpu_speedup_vs_n <- ggplot(subset(speedups,reconfig!=0), aes((n+1)*(n+1), solve)) +
                                    geom_smooth(colour=col_sprs, formula = y ~ log(x)) + geom_point(colour=col_sprs)

speedups = get_speedups(gpu_Tesla_dense_avgs, cpu_dense_avgs)
solve_dense_cpu_speedup_vs_n <- ggplot(subset(speedups,reconfig!=0), aes((n+1)*(n+1), solve)) + 
                                    geom_smooth(colour=col_dns, formula = y~log(x)) + geom_point(colour=col_dns)

speedups = get_speedups(gpu_Tesla_sparse_avgs, gpu_Tesla_dense_avgs)
solve_gpus_speedup_vs_n <- ggplot(subset(speedups,reconfig!=0 & n < 200), aes((n+1)*(n+1), solve)) + 
                                  geom_smooth(colour=col_dnsspr, formula = y~log(x)) + geom_point(colour=col_dnsspr)

## femses 

speedups = get_speedups(gpu_Tesla_femses_avgs, cpu_sparse_avgs)
solve_femses_sparse_cpu_speedup_vs_n <- ggplot(subset(speedups,reconfig!=0), aes((n+1)*(n+1), solve)) +
                              geom_smooth(colour=col_sprs, formula = y~log(x)) + geom_point(colour=col_sprs)
speedups = get_speedups(gpu_Tesla_femses_avgs, cpu_dense_avgs)
solve_femses_dense_cpu_speedup_vs_n <- ggplot(subset(speedups,reconfig!=0), aes((n+1)*(n+1), solve)) + 
                              geom_smooth(colour=col_dns, formula = y~log(x)) + geom_point(colour=col_dns)
speedups = get_speedups(gpu_Tesla_femses_avgs, gpu_Tesla_sparse_avgs)
solve_femses_sparse_gpu_speedup_vs_n <- ggplot(subset(speedups,reconfig!=0), aes((n+1)*(n+1), solve)) + 
                              geom_smooth(colour=col_dnsspr,formula = y~log(x)) + geom_point(colour=col_dnsspr)
speedups = get_speedups(gpu_Tesla_femses_avgs, gpu_Tesla_dense_avgs)
solve_femses_dense_gpu_speedup_vs_n <- ggplot(subset(speedups,reconfig!=0), aes((n+1)*(n+1), solve)) + 
                              geom_smooth(colour=col_femses, formula = y ~ log(x)) + geom_point(colour=col_femses)

## figures

solve_sparse_cpu_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                      labs(title="Speedup Over Serial Sparse Solver Time vs Problem Size", x ="Degrees of Freedom", y = "Speedup")
solve_dense_cpu_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                      labs(title="Speedup Over Serial Dense Solver Time vs Problem Size", x ="Degrees of Freedom", y = "Speedup")
solve_gpus_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                      labs(title="Speedup of GPU Sparse Solver Over GPU Dense Solver Time vs Problem Size", 
                        x ="Degrees of Freedom", y = "Speedup")
solve_femses_sparse_cpu_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                      labs(title="FEMSES Speedup Over Serial Sparse Solver Time vs Problem Size", x ="Degrees of Freedom", y = "Speedup")
solve_femses_dense_cpu_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                      labs(title="FEMSES Solver Speedup Over Serial Dense Solver Time vs Problem Size", x ="Degrees of Freedom", y = "Speedup")
solve_femses_sparse_gpu_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                      labs(title="FEMSES Solver Speedup Over GPU Sparse Solver Time vs Problem Size", x ="Degrees of Freedom", y = "Speedup")
solve_femses_dense_gpu_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                      labs(title="FEMSES Solver Speedup Over GPU Dense Solver Time vs Problem Size", x ="Degrees of Freedom", y = "Speedup")


###### Plot of allocation times #####
## add 4th plot with speedups of FEMSES & SPARSE

alloc_dense_vs_n = ggplot(subset(gpu_Tesla_dense_avgs,reconfig!=0), aes((n+1)*(n+1), allocation)) + 
                    geom_smooth(colour=col_sprs, formula = y ~ x*x) + geom_point(colour=col_sprs)
alloc_sparse_vs_n = ggplot(subset(gpu_Tesla_sparse_avgs,reconfig!=0), aes((n+1)*(n+1), allocation)) +
                    geom_smooth(colour=col_dns) + geom_point(colour=col_dns)
alloc_femses_vs_n = ggplot(subset(gpu_Tesla_femses_avgs,reconfig!=0), aes((n+1)*(n+1), allocation)) +
                    geom_smooth(colour=col_femses) + geom_point(colour=col_femses)


speedups = get_speedups(gpu_Tesla_sparse_avgs, gpu_Tesla_dense_avgs)
alloc_sparse_dense_speedup_vs_n = ggplot(subset(speedups,reconfig!=0 & n < 200), aes((n+1)*(n+1), allocation)) +
                    geom_smooth(colour=col_sprs, formula = y ~ log(x)) + geom_point(colour=col_sprs)
speedups = get_speedups(gpu_Tesla_femses_avgs, gpu_Tesla_sparse_avgs)
alloc_femses_sparse_speedup_vs_n = ggplot(subset(speedups,reconfig!=0 & n < 200), aes((n+1)*(n+1), allocation)) +
                    geom_smooth(colour=col_dns, formula = y ~ log(x)) + geom_point(colour=col_dns)
speedups = get_speedups(subset(gpu_Tesla_femses_avgs,allocation < 5), gpu_Tesla_dense_avgs)
alloc_femses_dense_speedup_vs_n = ggplot(subset(speedups,reconfig!=0 & n < 200), aes((n+1)*(n+1), allocation)) + 
                    geom_smooth(colour=col_femses, formula = y ~ log(x)) + geom_point(colour=col_femses)


alloc_dense_vs_n +  labs(title="Average Allocation Times for Dense Solver vs Problem Size", x ="Degrees of Freedom", y = "Allocation (ms)")
alloc_sparse_vs_n + labs(title="Average Allocation Times for Sparse Solver vs Problem Size", x ="Degrees of Freedom", y = "Allocation (ms)")
alloc_femses_vs_n + labs(title="Average Allocation Times for FEMSES Solver vs Problem Size", x ="Degrees of Freedom", y = "Allocation (ms)")
alloc_sparse_dense_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                      labs(title="GPU Sparse Speedup Over GPU Dense Allocation Time vs Problem Size", 
                        x ="Degrees of Freedom", y = "Speedup")
alloc_femses_sparse_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                      labs(title="FEMSES Speedup Over GPU Sparse Allocation Time vs Problem Size", x ="Degrees of Freedom", y = "Speedup")
alloc_femses_dense_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                      labs(title="FEMSES Speedup Over GPU Dense Allocation Time vs Problem Size", x ="Degrees of Freedom", y = "Speedup")


###### Plot of transfer times #####
## add 4th plot with speedups of FEMSES & SPARSE

transf_dense_vs_n = ggplot(subset(gpu_Tesla_dense_avgs,reconfig!=0), aes((n+1)*(n+1), transfer)) + 
                              geom_smooth(colour=col_dns) + geom_point(colour=col_dns)
transf_sparse_vs_n = ggplot(subset(gpu_Tesla_sparse_avgs,reconfig!=0), aes((n+1)*(n+1), transfer)) + 
                              geom_smooth(colour=col_sprs) + geom_point(colour=col_sprs)

speedups = get_speedups(gpu_Tesla_sparse_avgs, gpu_Tesla_dense_avgs)
transf_sparse_dense_speedup_vs_n = ggplot(subset(speedups,reconfig!=0 & n < 200), aes((n+1)*(n+1), transfer)) +
                              geom_smooth(colour=col_femses, formula = y ~ log(x)) + geom_point(colour=col_femses)


transf_dense_vs_n + labs(title="Average Transfer Times for Dense Solver vs Problem Size", x ="Degrees of Freedom", y = "Transfer (ms)")
transf_sparse_vs_n + labs(title="Average Transfer Times for Sparse Solver vs Problem Size", x ="Degrees of Freedom", y = "Transfer (ms)")
transf_sparse_dense_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                      labs(title="Sparse Speedup Over Dense Transfer Time vs Problem Size", x ="Degrees of Freedom", y = "Speedup")


###### Plots of element matrices evaluation times #####

speedups = get_speedups(gpu_Tesla_sparse_avgs, cpu_sparse_avgs)
elem_mats_dev_cpu_speedup_vs_b <- ggplot(subset(speedups,reconfig!=0 & n>=99), aes(block_size_X, elem_mats, 
                                                  colour = as.factor((n+1)*(n+1)))) + geom_smooth(formula = y~log(x)) + geom_point()
elem_mats_dev_cpu_speedup_vs_n <- ggplot(subset(speedups,reconfig!=0 & n>=99), aes((n+1)*(n+1), elem_mats)) +
                              geom_smooth(colour=col_femses, formula = y~log(x)) + geom_point(colour=col_femses)

labels = unique(subset(speedups,n>=99)[1])
labels = list((labels+1) **2 )

n_labeller <- function(variable,value){
  return(labels[value])
}


elem_mats_dev_speedups_reconfig <- ggplot(subset(speedups,n>=99 & block_size_X < 224), aes(block_size_X, elem_mats,
                                        colour = as.factor(reconfig))) + geom_point() + geom_smooth(formula = y ~ log(x)) + 
                                        facet_wrap(.~ n, scales="free", labeller = as_labeller(labels))


elem_mats_dev_cpu_speedup_vs_b + scale_y_log10(breaks = base_breaks()) + 
                                labs(title="Speedup of GPU Element Matrices Generation Over CPU vs Block Size", x ="Block Size", y = "Speedup") +
                                scale_color_discrete(name="DOF")
elem_mats_dev_cpu_speedup_vs_n + scale_y_log10(breaks = base_breaks()) + 
                                labs(title="Speedup of GPU Element Matrices Generation Over CPU vs Problem Size",x ="Degrees of Freedom", y = "Speedup")

elem_mats_dev_speedups_reconfig + scale_y_log10(breaks = base_breaks()) +
                                labs(title="Speedup of GPU Element Matrices Generation Over CPU vs Block Size for\nIndividual Problem Sizes", 
                                x ="Block Size", y = "Speedup") +
                                scale_color_discrete(name="Mem Config", labels =c("1:3","1:1","3:1"))



###### Plots of assembly times #####

## sparse
speedups = get_speedups(gpu_Tesla_sparse_avgs, cpu_sparse_avgs)
assem_dev_cpu_sparse_speedup_vs_b <- ggplot(subset(speedups,reconfig!=0 & n>=99), aes(block_size_X, assembly,
                                            colour = as.factor((n+1)*(n+1)))) + geom_point() + geom_smooth(formula = y ~ log(x))
assem_dev_cpu_sparse_speedup_vs_n <- ggplot(subset(speedups,reconfig!=0 & n>99), aes((n+1)*(n+1), assembly)) +
                              geom_smooth(colour=col_dns, formula = y~log(x)) + geom_point(colour=col_dns)

assem_dev_sparse_speedups_reconfig <- ggplot(subset(speedups,n>=99 & block_size_X < 224), aes(block_size_X, assembly,
                                          colour = as.factor(reconfig))) + geom_point() + geom_smooth(formula = y ~ log(x)) +
                                          facet_wrap(~ n, scales="free")

## dense
speedups = get_speedups(gpu_Tesla_dense_avgs, cpu_dense_avgs)
assem_dev_cpu_dense_speedup_vs_b <- ggplot(subset(speedups,reconfig!=0 & n>10), aes(block_size_X, assembly,
                                             colour = as.factor((n+1)*(n+1)))) + geom_point() + geom_smooth(formula = y ~ log(x))
assem_dev_cpu_dense_speedup_vs_n <- ggplot(subset(speedups,reconfig!=0), aes((n+1)*(n+1), assembly)) +
                              geom_smooth(colour=col_sprs,formula = y ~ log(x)) + geom_point(colour=col_sprs)

assem_dev_dense_speedups_reconfig <- ggplot(subset(speedups,n>10 & block_size_X < 224), aes(block_size_X, assembly,
                                        colour = as.factor(reconfig))) + geom_point() + geom_smooth(formula = y ~ log(x)) +
                                        facet_wrap(~ n, scales="free")

## sparse vs dense
speedups = get_speedups(gpu_Tesla_sparse_avgs, gpu_Tesla_dense_avgs)
assem_dev_gpu_sparse_dense_speedup_vs_b <- ggplot(subset(speedups,reconfig!=0 & n <= 199), aes(block_size_X, assembly,
                                 colour = as.factor((n+1)*(n+1)))) + geom_point() + geom_smooth(formula = y ~ log(x))


## figures
assem_dev_cpu_sparse_speedup_vs_b + scale_y_log10(breaks = base_breaks()) + 
                                labs(title="Speedup of GPU Sparse Stiffness Matrix Assembly Over CPU vs Block Size", x ="Block Size", y = "Speedup") +
                                scale_color_discrete(name="DOF")
assem_dev_cpu_sparse_speedup_vs_n + scale_y_log10(breaks = base_breaks()) + 
                                labs(title="Speedup of GPU Sparse Stiffness Matrix Assemble Over CPU vs Problem Size",
                                  x ="Degrees of Freedom", y = "Speedup")
assem_dev_sparse_speedups_reconfig + scale_y_log10(breaks = base_breaks()) + 
                                labs(title="Speedup of GPU Sparse Stiffness Matrix Assembly Over CPU vs Block Size for\nIndividual Problem Sizes", 
                                x ="Block Size", y = "Speedup") +
                                scale_color_discrete(name="Mem Config", labels =c("1:3","1:1","3:1"))

assem_dev_cpu_dense_speedup_vs_b + scale_y_log10(breaks = base_breaks()) + 
                                labs(title="Speedup of GPU Dense Stiffness Matrix Assembly Over CPU vs Block Size", x ="Block Size", y = "Speedup") +
                                scale_color_discrete(name="DOF")
assem_dev_cpu_dense_speedup_vs_n + labs(title="Speedup of GPU Dense Stiffness Matrix Assemble Over CPU vs Problem Size",
                                  x ="Degrees of Freedom", y = "Speedup")
assem_dev_dense_speedups_reconfig + scale_y_log10(breaks = base_breaks()) + 
                                labs(title="Speedup of GPU Dense Stiffness Matrix Assembly Over CPU vs Block Size for\nIndividual Problem Sizes", 
                                x ="Block Size", y = "Speedup") +
                                scale_color_discrete(name="Mem Config", labels =c("1:3","1:1","3:1"))

assem_dev_gpu_sparse_dense_speedup_vs_b + scale_y_log10(breaks = base_breaks()) + 
                                labs(title="Speedup of GPU Sparse Stiffness Matrix Assembly Over GPU Dense\nAssembly vs Block Size", 
                                x ="Block Size", y = "Speedup") +
                                scale_color_discrete(name="DOF")
# trade off between solver time and assembly time #


###### Plots of element matrices + assembly  #####

## sparse
speedups = get_speedups(gpu_Tesla_sparse_avgs, cpu_sparse_avgs)
elems_p_assem_ker_cpu_sparse_speedup_vs_b <- ggplot(subset(speedups, reconfig != 0 & n > 50), aes(block_size_X, elems_p_assemb,
                                              colour = as.factor((n+1)*(n+1)))) + geom_point() + geom_smooth(formula = y ~ log(x))
elems_p_assem_ker_cpu_sparse_speedup_vs_n <- ggplot(subset(speedups, reconfig != 0 & n > 50), aes((n+1)*(n+1), elems_p_assemb)) + 
                                              geom_point(colour = col_dnsspr) + geom_smooth(formula = y ~ log(x), colour = col_dnsspr)

elems_p_assem_ker_sparse_speedups_reconfig <- ggplot(subset(speedups,block_size_X < 224 & n > 50), aes(block_size_X, elems_p_assemb,
                                              colour = as.factor(reconfig))) + geom_point() + geom_smooth(formula = y~log(x)) +
                                              facet_wrap(~ n, scales="free")

## dense
speedups = get_speedups(gpu_Tesla_dense_avgs, cpu_dense_avgs)
elems_p_assem_ker_cpu_dense_speedup_vs_b <- ggplot(subset(speedups,reconfig!=0 & n > 10), aes(block_size_X, elems_p_assemb,
                                                                  colour = as.factor((n+1)*(n+1)))) + geom_point() +
                                                                  geom_smooth(formula = y~log(x))
elems_p_assem_ker_cpu_dense_speedup_vs_n <- ggplot(subset(speedups, reconfig != 0 & n > 10), aes((n+1)*(n+1), elems_p_assemb)) + 
                                              geom_point(colour = col_femses) + geom_smooth(colour = col_femses)

elems_p_assem_ker_dense_speedups_reconfig <- ggplot(subset(speedups,block_size_X < 224 & n > 10), aes(block_size_X, elems_p_assemb,
                                              colour = as.factor(reconfig))) + geom_point() + geom_smooth(formula = y~log(x)) +
                                              facet_wrap(~ n, scales="free") 

## sparse vs dense
speedups = get_speedups(gpu_Tesla_sparse_avgs, gpu_Tesla_dense_avgs)
elems_p_assem_ker_gpu_sparse_dense_speedup_vs_b <- ggplot(subset(speedups,reconfig!=0 & n <= 199), aes(block_size_X, elems_p_assemb,
                                                                  colour = as.factor((n+1)*(n+1)))) + geom_point() +
                                                                  geom_smooth(formula = y~log(x))


## figures
elems_p_assem_ker_cpu_sparse_speedup_vs_b + scale_y_log10(breaks = base_breaks()) +
                                labs(title="Speedup of GPU Sparse Main Kernel Over CPU vs Block Size", x ="Block Size", y = "Speedup") +
                                scale_color_discrete(name="DOF")
elems_p_assem_ker_cpu_sparse_speedup_vs_n + scale_y_log10(breaks = base_breaks()) +
                                labs(title="Speedup of GPU Sparse Main Kernel Over CPU vs Problem Size", x ="Degrees of Freedom", y = "Speedup")

elems_p_assem_ker_sparse_speedups_reconfig + scale_y_log10(breaks = base_breaks()) + 
                                labs(title="Speedup of GPU Sparse Main Kernel Over CPU vs Block Size for\n Individual Problem Sizes", 
                                x ="Block Size", y = "Speedup") +
                                scale_color_discrete(name="Mem Config", labels =c("1:3","1:1","3:1")) 

elems_p_assem_ker_cpu_dense_speedup_vs_b + scale_y_log10(breaks = base_breaks()) + 
                                labs(title="Speedup of GPU Dense Main Kernel Over CPU vs Block Size", 
                                x ="Block Size", y = "Speedup") +
                                scale_color_discrete(name="DOF")
elems_p_assem_ker_cpu_dense_speedup_vs_n +  labs(title="Speedup of GPU Dense Main Kernel Over CPU vs Problem Size",
                                                x ="Degrees of Freedom", y = "Speedup")
elems_p_assem_ker_dense_speedups_reconfig + scale_y_log10(breaks = base_breaks()) + 
                                labs(title="Speedup of GPU Dense Main Kernel Over CPU vs Block Size for\n Individual Problem Sizes", 
                                x ="Block Size", y = "Speedup") +
                                scale_color_discrete(name="Mem Config", labels =c("1:3","1:1","3:1")) 

elems_p_assem_ker_gpu_sparse_dense_speedup_vs_b + scale_y_log10(breaks = base_breaks()) + 
                                labs(title="Speedup of GPU Sparse Kernel Over Dense Kernel vs Block Size", x ="Block Size", y = "Speedup") +
                                scale_color_discrete(name="DOF")


###### Conversion time, sparsity scan & num_iterations vs problem size #####

conv_vs_n = ggplot(gpu_Tesla_dnsspr_avgs, aes((n+1)*(n+1), convert)) + geom_smooth(colour=col_sprs) + geom_point(colour=col_sprs)
sparsity_v_n = ggplot(cpu_sparse_avgs, aes((n+1)*(n+1), sparsity.scan)) + geom_smooth(colour=col_dns) + geom_point(colour=col_dns)
iters_v_n = ggplot(subset(gpu_Tesla_femses_avgs,block_size_X==32 & reconfig==0), aes((n+1)*(n+1), iterations)) + 
                    geom_segment(aes(x = 0, y = 0, xend = 40000, yend = 40000),linetype="dashed", colour=col_sprs, size=1.0) + 
                    geom_smooth(colour=col_femses) + geom_point(colour=col_femses)

conv_vs_n + labs(title="Conversion from Dense to CSR Time Taken by GPU vs. Problem Size", x ="Degrees of Freedom", y = "Conversion (ms)")
sparsity_v_n + labs(title="Time Taken to Perform Sparsity Pass by CPU vs Problem Size", x ="Degrees of Freedom", y = "Sparsity Pass (ms)")
iters_v_n + labs(title="Number of Iterations for FEMSES to Converge vs. Problem Size", x ="Degrees of Freedom", y = "Iterations")



###### CPU Profiling ######

total_sparse_cpu = ggplot(cpu_sparse_avgs, aes((n+1)*(n+1), total)) + 
                    geom_smooth(colour=col_sprs) + geom_point(colour=col_sprs)
total_dense_cpu = ggplot(cpu_dense_avgs, aes((n+1)*(n+1), total)) + 
                    geom_smooth(colour=col_dns) + geom_point(colour=col_dns)
speedups = get_speedups(cpu_sparse_avgs, cpu_dense_avgs)
total_sparse_dense_cpu_speedup_vs_n <- ggplot(speedups, aes((n+1)*(n+1), total)) + 
                              geom_smooth(colour=col_femses, formula = y ~ log(x)) + geom_point(colour=col_femses)

elem_mats_cpu = ggplot(rbind(cpu_sparse_avgs, cpu_dense_avgs), aes((n+1)*(n+1), elem_mats)) + 
                    geom_smooth(colour=col_dns) + geom_point(colour=col_dns)
assembly_sparse_cpu = ggplot(cpu_sparse_avgs, aes((n+1)*(n+1), assembly)) + 
                    geom_smooth(colour=col_femses) + geom_point(colour=col_femses)
assembly_dense_cpu = ggplot(cpu_dense_avgs, aes((n+1)*(n+1), assembly)) + 
                    geom_smooth(colour=col_dnsspr) + geom_point(colour=col_dnsspr)
solve_sparse_cpu = ggplot(cpu_sparse_avgs, aes((n+1)*(n+1), solve)) + 
                    geom_smooth(colour=col_sprs) + geom_point(colour=col_sprs)
solve_dense_cpu = ggplot(cpu_dense_avgs, aes((n+1)*(n+1), solve)) + 
                    geom_smooth(colour=col_dns) + geom_point(colour=col_dns)
sparsity_pass_cpu = ggplot(cpu_sparse_avgs, aes((n+1)*(n+1), sparsity.scan)) + 
                    geom_smooth(colour=col_sprs) + geom_point(colour=col_sprs)

data.m <- melt(subset(cpu_sparse_avgs, n > 25), id.vars=c("n", "m", "total", "block_size_X", "reconfig", "sse", 
                                "iterations", "elems_p_assemb", "convert", "transfer"))
prop_sparse <- ggplot(data.m, aes(fill=variable, y=value/total, x=as.factor((n+1)*(n+1)))) +
                                geom_bar(position="fill", stat="identity") +
                                labs(title="Proportion of Computation Time Taken By Each Step\nin Sparse Serial Process", 
                                x ="Problem Size", y = NULL, fill = "Kernel") +
                                scale_y_continuous(labels = scales::percent_format()) + 
                                scale_fill_brewer(palette="Set2")

data.m <- melt(subset(cpu_dense_avgs, n > 10), id.vars=c("n", "m", "total", "block_size_X", "reconfig", "sse", 
                                "iterations", "elems_p_assemb", "convert", "transfer", "sparsity.scan"))
prop_dense <- ggplot(data.m, aes(fill=variable, y=value/total, x=as.factor((n+1)*(n+1)))) +
                                geom_bar(position="fill", stat="identity") +
                                labs(title="Proportion of Computation Time Taken By Each Step\nin Dense Serial Process", 
                                x ="Problem Size", y = NULL, fill = "Kernel") + 
                                scale_y_continuous(labels = scales::percent_format()) + 
                                scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest2"))


total_sparse_cpu +  labs(title="Total Computation Time for Sparse CPU Solution vs Problem Size",
                          x ="Degrees of Freedom", y = "Total (ms)")
total_dense_cpu +  labs(title="Total Computation Time for Dense CPU Solution vs Problem Size", 
                          x ="Degrees of Freedom", y = "Total (ms)")
total_sparse_dense_cpu_speedup_vs_n + labs(title="Speedup of Total Time for Sparse Solution\nover Dense vs Problem Size", 
                          x ="Degrees of Freedom", y = "Speedup")

solve_sparse_cpu +  labs(title="Sparse Solver Times for CPU vs Problem Size", 
                          x ="Degrees of Freedom", y = "Solve (ms)")
solve_dense_cpu +  labs(title="Dense Times for CPU vs Problem Size", 
                          x ="Degrees of Freedom", y = "Solve (ms)")
elem_mats_cpu +  labs(title="Element Matrices Generation Times for CPU vs Problem Size", 
                          x ="Degrees of Freedom", y = "Element Matrices Generation (ms)")
assembly_sparse_cpu +  labs(title="Sparse Global Stiffness Matrix Assembly Times for CPU\nvs Problem Size", 
                          x ="Degrees of Freedom", y = "Assembly (ms)")
assembly_dense_cpu +  labs(title="Dense Global Stiffness Matrix Assembly Times for CPU\nvs Problem Size", 
                          x ="Degrees of Freedom", y = "Assembly (ms)")
sparsity_pass_cpu +  labs(title="Sparsity Pass Times vs Problem Size", 
                          x ="Degrees of Freedom", y = "Sparsity Pass (ms)")

prop_sparse
prop_dense




