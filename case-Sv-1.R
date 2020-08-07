library(rstan)
library(ggmcmc)

#======  input ===============================#
input_file_name <- 'case-Sv-1'
pars <- c('s11', 's12', 's13', 's22', 's23', 's33', 'p0', 'K', 'G', 'tv', 'ts', 'sigma', 'eigen_values', 'eigen_vectors')

#compile the model?
compile <- "yes"

#parallelize chains ?
parallel <- "yes"

# directional cosine vectors
xx1 <- 1.0
xx2 <- 0.0
xx3 <- 0.0

yy1 <- 0.0
yy2 <- 1.0
yy3 <- 0.0

zz1 <- 0.0
zz2 <- 0.0
zz3 <- 1.0

xy1 <- cos(pi/4.0)
xy2 <- sin(pi/4.0)
xy3 <- 0.0

yx1 <- cos(pi/4.0)
yx2 <- -sin(pi/4.0)
yx3 <- 0.0

yz1 <- 0.0
yz2 <- cos(pi/4.0)
yz3 <- sin(pi/4.0)

zy1 <- 0.0
zy2 <- cos(pi/4.0)
zy3 <- -sin(pi/4.0)

zx1 <- sin(pi/4.0)
zx2 <- 0.0
zx3 <- cos(pi/4.0)

xz1 <- -sin(pi/4.0)
xz2 <- 0.0
xz3 <- cos(pi/4.0)

n <- list(c(xx1, xx2, xx3), 
          c(yy1, yy2, yy3), 
          c(zz1, zz2, zz3), 
          c(xy1, xy2, xy3), 
          c(yx1, yx2, yx3), 
          c(yz1, yz2, yz3),
          c(zy1, zy2, zy3),
          c(zx1, zx2, zx3),
          c(xz1, xz2, xz3)
          )
#================================================#


#==== don't change the below section!!===========================#
input_file <- paste("input/", input_file_name, ".csv",sep="")

output <- paste("output/", input_file_name, '.RData', sep="")
summary_name <- paste("output/", input_file_name, '_summary.txt', sep="")
plot_name <- paste("output/", input_file_name, '.pdf', sep="")

input <-read.csv(input_file, stringsAsFactors = F)
seed <- input$seed
chains <- input$chains
iter <- input$iteration
warmup <- input$warmup
thin <- input$thin
cores <- switch(parallel,
                "yes" = chains,
                "no" = 1L,
                stop("please input yes or no for parallelization!")
                )
used_CH <- c(input$xx, 
             input$yy, 
             input$zz, 
             input$xy, 
             input$yx, 
             input$yz, 
             input$zy, 
             input$zx, 
             input$xz
             )

Sv <- input$Sv
dt <- input$delay_time.h.


#chNo <- select(input, num_range(prefix = "ch", range = 1:channel_No,1))

data_name <- paste("data/", input$data.file, sep="")
model <- paste("model/", input$model, sep="")
model_name <- paste(model, ".stan", sep="") 
compiled_model_name <- paste(model, ".rds", sep="")

d <- read.csv(file = data_name)
N <- nrow(d)
D <- length(used_CH)
I <- diag(D)

#===== data used for modeling =========#
t <- d$time
en <- d[, used_CH]
#======================================#

data <- list(N = N, D = D, t = t, n = n, I = I, en = en, dt = dt, Sv = Sv)


rstan_options(auto_write = TRUE)
stanmodel <- switch (compile,
                   "yes" = stan_model(file = model_name),
                   "no" =  readRDS(compiled_model_name),
                   stop("please input yes or no for compile!")
                   )


fit <- sampling(
  stanmodel,
  data = data,
  pars = pars,
  seed = seed,
  chains = chains,
  cores = cores,
  iter = iter,
  warmup = warmup,
  thin = thin
  #, init=function(){
  #  list(K=runif(1,-10,10),G=runif(1,0,10))
  #}
)

save.image(file = output)

print("================ plot information ===============")

write.table(data.frame(summary(fit)$summary),
            file = summary_name,
            sep = '\t',
            quote = FALSE,
            col.names = NA
            )

 
pdf(file = plot_name)

ggmcmc(ggs(fit, inc_warmup = T, stan_include_auxiliar = T),
       file = NULL, 
       width = 7,
       height = 7,
       param_page = 4,
       plot = "traceplot"
       )

ggmcmc(ggs(fit, stan_include_auxiliar = T),
       file = NULL, 
       width = 7,
       height = 7,
       param_page = 3
       )

ggs_pairs(ggs(fit), lower = list(continuous = "density"))

dev.off()
