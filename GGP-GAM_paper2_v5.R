#' title: "Encapsulating spatially varying relationships with a Generalized Additive Model"

#' # Abstract
#' 
#' # 1. Introduction 
#' 
#' # 2. Literature Review
#' ## 2.1 SVC models
#' ## GAM based SVC models
#' 
#' # 3. Methodology
#' ## 3.1 Data
#' 
# Use DMs ESF matrices data to creategrids
library(tidyverse)
library(spmoran)
library(sf)
gnum<-25
px    <-rep(1:gnum,gnum)
py    <-rep(1:gnum,each=gnum)
coords<-cbind(px,py)
meig<-meigen(coords)

b0 = 2
b1 <- meig$sf[,1] %>% scale()
b2<-meig$sf[,10] %>% scale()
b3<-meig$sf[,25] %>% scale()
b123 = b1+b2+b3
gr <-data.frame(coords,b0, b1,b2,b3)

# Figure 1
library(tidyverse)
library(cowplot)
library(cols4all)
# c4a_gui() # uncomment to see cols4all() palettes
# function to make the maps with ggplot 
beta_plot_function = function(gr, beta, tit, legend.pos = "none", 
                              lims = c(-2.2, 2.2), c4apal = "kovesi.bu_yl_mg") {
  ggplot(gr, aes(x = px, y = py, fill = beta)) +
  geom_raster()+
  coord_fixed() +
  scale_fill_continuous_c4a_seq(palette=c4apal, limits = lims, reverse = T) + 
  theme_bw() + xlab("") + ylab("") +  
  ggtitle(tit) +
  theme(legend.position = legend.pos, 
        legend.title = element_blank(),
        legend.key.width = unit(.5, "cm"),
        legend.key.height = unit(.65, "cm"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 
}
# create the plots 
tit =expression(paste("True "*beta[0]*""))
p0 = beta_plot_function(gr, beta = gr$b0, tit) 
tit =expression(paste("True "*beta[1]*""))
p1 = beta_plot_function(gr, beta = gr$b1, tit) 
tit =expression(paste("True "*beta[2]*""))
p2 = beta_plot_function(gr, beta = gr$b2, tit) 
tit =expression(paste("True "*beta[3]*""))
p3 = beta_plot_function(gr, beta = gr$b3, tit,legend.pos = "right") 
legend_betas <- get_legend(p3)
p3 = p3 + theme(legend.position = "none")
# and now plot together using plot_grid from cowplot
plot_grid(p0, p1, p2, p3, legend_betas, nrow = 1, rel_widths = c(1,1,1,1,0.25))

## Make simulated data
# Loads precompiled results 
# Edit MakeSim = T if you want to create from fresh 
# The individual simulated datastes will be different 
MakeSim = F
if(MakeSim) {
  ## Create Simulations
  make_x = function(){
  	x = rnorm(nrow(gr)) 
  	x = (x - min(x))/(max(x) - min(x))
  	x
  }	
  y_sim = x1_sim = x2_sim = x3_sim = e_sim = vector()
  for(i in 1:100){
  	x1 = make_x()
  	x2 = make_x()
  	x3 = make_x()
  	e = make_x()*.25
  	y.i = gr$b0 + (gr$b1*x1) + (gr$b2*x2) + (gr$b3*x3) + e
  	y_sim = cbind(y_sim, y.i)
  	x1_sim= cbind(x1_sim, x1)
  	x2_sim= cbind(x2_sim, x2)
  	x3_sim= cbind(x3_sim, x3)
  	e_sim= cbind(e_sim, e)
  }
  save(list = c("y_sim", "x1_sim", "x2_sim", "x3_sim", "e_sim"), file = "sim_paper2.RData")
} else {
  load("sim_paper2.RData")
}

# Load packages and define function 
library(tidyverse)
library(spmoran)
library(sf)
library(mgcv)
library(parallel)  
library(doParallel)
# function to exract ESF matrices data to create grids
make_esf = function(gnum, big = T) {
	px    <-rep(1:gnum,gnum)
	py    <-rep(1:gnum,each=gnum)
	coords<-cbind(px,py)
	if(!big) meig<-meigen(coords) 
	if(big)  meig<-meigen_f(coords) # for large samples (e.g., n > 5000)
	b0 = 2
	b1 <- meig$sf[,1] %>% scale()
	b2<-meig$sf[,10] %>% scale()
	b3<-meig$sf[,25] %>% scale()
	b123 = b1+b2+b3	
	gr <-data.frame(coords,b0, b1,b2,b3)
	# create simulations
	make_x = function(){
	  	x = rnorm(nrow(gr)) 
	  	x = (x - min(x))/(max(x) - min(x))
	  	x
	}	
	y_sim = x1_sim = x2_sim = x3_sim = e_sim = vector()
	#for(i in 1:100){
	for(i in 50){				# a single instance
	  	x1 = make_x()
	  	x2 = make_x()
	  	x3 = make_x()
	  	e = make_x()*.25
	  	y.i = gr$b0 + (gr$b1*x1) + (gr$b2*x2) + (gr$b3*x3) + e
	  	y_sim = cbind(y_sim, y.i)
	  	x1_sim= cbind(x1_sim, x1)
	  	x2_sim= cbind(x2_sim, x2)
	  	x3_sim= cbind(x3_sim, x3)
	}
	return(list(y_sim = y_sim, x1_sim = x1_sim, x2_sim = x2_sim, x3_sim = x3_sim,
				px = px, py = py, gr = gr))
}
# Create data and save (use of rnorm in the above function)
# Loads precompiled results 
# Edit MakeBigSim = T if you want to create from fresh 
MakeBigSim = F
if (MakeBigSim) {
  gnum = 500
  sims = make_esf(gnum, big = T)
  save(sims, file = "sims500.RData")
} else {
  load("sims500.RData")
}

# Figure 2
# define gr from loaded data
gr_big = sims$gr
lims = c(-3, 3)
tit =expression(paste("True "*beta[0]*""))
p0b = beta_plot_function(gr_big, beta = gr_big$b0, tit, lims = lims,
                        c4apal = "kovesi.bu_yl_mg") 
tit =expression(paste("True "*beta[1]*""))
p1b = beta_plot_function(gr_big, beta = gr_big$b1, tit, lims = lims,
                        c4apal = "kovesi.bu_yl_mg")  
tit =expression(paste("True "*beta[2]*""))
p2b = beta_plot_function(gr_big, beta = gr_big$b2, tit, lims = lims, 
                        c4apal = "kovesi.bu_yl_mg")  
tit =expression(paste("True "*beta[3]*""))
p3b = beta_plot_function(gr_big, beta = gr_big$b3, tit, lims = lims,
                        legend.pos = "right",
                        c4apal = "kovesi.bu_yl_mg")  
legend_betas <- get_legend(p3)
p3b = p3b + theme(legend.position = "none")
# and now plot together using plot_grid from cowplot
plot_grid(p0b, p1b, p2b, p3b, legend_betas, nrow = 1, rel_widths = c(1,1,1,1,0.25))

#' ## 3.2 GAM-based SVC models
#' ## 3.3 Analysis I: comparing GGP-GAM and MGWR

# Figure 3
library(mgcv)
library(tidyverse)
# from Clark 2017
# see https://m-clark.github.io/generalized-additive-models/technical.html#the-number-of-knots-and-where-to-put-them
set.seed(123)
x  = runif(500)
mu = sin(2 * (4 * x - 2)) + 2 * exp(-(16 ^ 2) * ((x - .5) ^ 2))
y  = rnorm(500, mu, .3)
res = NULL
for(i in c(3, 6, 9, 15, 21)){
	gam.i = gam(y~s(x, bs = 'gp', k = i))
	y.i = gam.i$fitted.values
	res = cbind(res, y.i)
}
res = cbind(x = x, res)
res = as_tibble(res)
colnames(res) = c("x", paste0("k=", c(3, 6, 9, 15, 21)))
res <- res |> pivot_longer(-x) |> 
	mutate(Knots = factor(name, levels = paste0("k=", c(3, 6, 9, 15, 21)))) 
p1a = ggplot() + 
	geom_point(data = data.frame(x,y), aes(x,y), col = "lightgrey", alpha = 0.7) +
	geom_line(data = res, aes(x, y = value, col = Knots), lwd = 1) +
  scale_color_discrete_c4a_cat(palette="tol.muted") +
	ylab("y") +
	theme_bw()
p1a

#' ## 3.4 Analysis II: GGP-GAM tuning with a larger dataset
#' ## 3.5 Empirical case study
#' 
# Figure 4
gb = st_read("data_gb.geojson", quiet = TRUE)
gb <- gb |> mutate(one=1) |> st_transform(27700)
gb_df <- st_centroid(gb) |> 
  st_coordinates() |> 
  as_tibble() |> 
  bind_cols(gb)

ggplot(gb,aes(fill=100*(share_leave - 0.5))) + 
  geom_sf(lwd = 0) + 
  ggthemes::theme_map() + 
  scale_fill_continuous_c4a_div(palette="tableau.red_blue_diverging",
                                name='Leave\n(%) - 50%' ) +
  theme(legend.position="right") + theme_map()

#' 
#' # 4. Results
#' ## 4.1 Comparing GGP-GAM and MGWR
#' 
library(mgcv)
library(GWmodel)
library(tidyverse)
library(spmoran)
library(sf)
library(parallel)  
library(doParallel)

## Functions - declare some functions
# r2
rsq <- function (x, y) cor(x, y) ^ 2
# rmse - how far the data from the line of best fit
rmse <- function(x, y) sqrt(mean((x - y)^2))
# mae = average absolute difference between obs and pred
mae <- function(x, y){
  n = length(x) 
	sum = 0
	for (i in 1:n){
		sum = abs(x[i] - y[i]) + sum
	}
	sum/n
}
# SVC extraction
ggp_gam_svc = function(model, terms, input_data) {
	n_t = length(terms)
	input_data_copy = input_data
	output_data = input_data
	for (i in 1:n_t) {
		zeros = rep(0, n_t)
		zeros[i] = 1
		df.j = vector()
		terms_df = data.frame(matrix(rep(zeros, nrow(input_data)), ncol = n_t, byrow = T))
		names(terms_df) = terms
		input_data_copy[, terms] = terms_df
		se.j = predict(model, se = TRUE, newdata = input_data_copy)$se.fit
		b.j=predict(model,newdata=input_data_copy)
		expr1 = paste0("b_", terms[i], "= b.j")
		expr2 = paste0("se_",terms[i], "= se.j")
		output_data = output_data %>% 
			mutate(within(., !!parse(text = expr1))) %>% 
			mutate(within(., !!parse(text = expr2))) 
	}
	output_data$yhat = predict(model, newdata = input_data)
	output_data
}

## GAM functions
# do the GAM
dogam_function = function(i){
	gam_res_tab = vector()
 	gb0.all = gb1.all = gb2.all = gb3.all = vector()
  	sps = ypred_gam  = vector()
  	terms = c("Intercept", "x1", "x2", "x3")
  	df.i = data.frame(y = y_sim[,i],
                      x1 = x1_sim[,i],
                      x2 = x2_sim[,i],
                      x3 = x3_sim[,i],
                      #e = e_sim[,i],
                      Intercept = 1, 
                      X = gr$px, Y = gr$py)
   	kval = 155
   	gam.i <- mgcv::gam(y~ 0 + 
                Intercept + s(X,Y,bs='gp',k=kval,by=Intercept) + 
                x1 + s(X,Y,bs='gp',k=kval,by=x1) +
                x2 + s(X,Y,bs='gp',k=kval,by=x2) + 
                x3 + s(X,Y,bs='gp',k=kval,by=x3), 
                data = df.i, 
                method = "REML")
    return(gam.i)
}

# extract the results	
gam_mods_accuracies = function(gam_mod){	
	gam.i = gam_mod
	df.i = data.frame(y = y_sim[,i],
                      x1 = x1_sim[,i],
                      x2 = x2_sim[,i],
                      x3 = x3_sim[,i],
                      #e = e_sim[,i],
                      Intercept = 1, 
                      X = gr$px, Y = gr$py)
    terms = c("Intercept", "x1", "x2", "x3")
    df.ii = ggp_gam_svc(gam.i, terms, df.i)
    ## evaluate betas prediction accuracy
    resids = gam.i$residuals
    aicval = nrow(df.i) * (log(2*pi)+1 + log( (sum(resids^2)/nrow(df.i) ) ) ) + ((4+1)*2)
    rsq.i = c(aicval, rsq(gr$b1, df.ii$b_x1), rsq(gr$b2, df.ii$b_x2), rsq(gr$b3, df.ii$b_x3))
    rmse.i = c(rmse(gr$b0, df.ii$b_Intercept), rmse(gr$b1, df.ii$b_x1), rmse(gr$b2, df.ii$b_x2), 
               rmse(gr$b3, df.ii$b_x3))
    mae.i = c(mae(gr$b0, df.ii$b_Intercept), mae(gr$b1, df.ii$b_x1), mae(gr$b2, df.ii$b_x2), 
              mae(gr$b3, df.ii$b_x3))
    gam_res = c(rsq.i, rmse.i, as.vector(mae.i))
    names(gam_res) = c("AIC", "Rsq_B1", "Rsq_B2", "Rsq_B3", "RMSE_B0", "RMSE_B1",
    		"RMSE_B2","RMSE_B3", "MAE_B0", "MAE_B1", "MAE_B2","MAE_B3")	
    ## model prediction accuracy
    gam.pred = as.vector(predict(gam.i, newdata = df.i))
    rsq.gam.y = rsq(gam.pred, df.i$y)
    rmse.gam.y = rmse(gam.pred, df.i$y)
    mae.gam.y = mae(gam.pred, df.i$y)
    ypred.gam = c(rsq.gam.y, rmse.gam.y, mae.gam.y)
    names(ypred.gam) = c("R2", "RMSE", "MAE")
    sps = as.vector(gam.i$sp)
    names(sps) = c("sp_gam_b0", "sp_gam_b1", "sp_gam_b2", "sp_gam_b3")
    return(list(gam_res = gam_res, gam_resids = resids, ypred_gam = ypred.gam, sps = sps))
}

## MGWR functions
# do the MGWR
dogwr_function = function(i){
    mgwr_res_tab = vector()
    mb0.all = mb1.all = mb2.all = mb3.all = vector()
    bws = ypred_mgwr = vector()
    # set up the data
    df.i = data.frame(y = y_sim[,i],
                      x1 = x1_sim[,i],
                      x2 = x2_sim[,i],
                      x3 = x3_sim[,i],
                      #e = e_sim[,i],
                      Intercept = 1, 
                      X = gr$px, Y = gr$py)
    df.i.sp = sp::SpatialPointsDataFrame(coords = df.i[, c("X", "Y")], 
                                     data = data.frame(df.i[, c("y", "x1", "x2", "x3")]))
	# find bandwidths
    bws.i = GWmodel::gwr.multiscale(y~x1+x2+x3,  
  	                        data = df.i.sp,
  	                        adaptive = T,
  	                        criterion="CVR",
  	                        kernel = "bisquare",
  	                        bws0=c(10,10,10,10),
  	                        verbose = F, predictor.centered=rep(T, 4))
    # now do with bandwidths to confirm
    msgwr_bws.i = as.vector(bws.i$GW.arguments$bws)
    bws.i = GWmodel::gwr.multiscale(y~x1+x2+x3,  
  	                        data = df.i.sp,
  	                        adaptive = T, 
  	                        criterion="CVR",
  	                        kernel = "bisquare",
                            bws0 = msgwr_bws.i,
                            bw.seled=rep(T, 4),
  	                        verbose = F, predictor.centered=rep(F, 4))
    return(bws.i)
}  
# extract the results
mgwr_mods_accuracies = function(mgwr_mod){	
	bws.i = mgwr_mod #mgwr_list[[1]]
	df.i = data.frame(y = y_sim[,i],
                      x1 = x1_sim[,i],
                      x2 = x2_sim[,i],
                      x3 = x3_sim[,i],
                      #e = e_sim[,i],
                      Intercept = 1, 
                      X = gr$px, Y = gr$py)
	## evaluate betas prediction accuracy
    resids = bws.i$SDF$residual
    aicval = nrow(df.i) * (log(2*pi)+1 + log( (sum(resids^2)/nrow(df.i) ) ) ) + ((4+1)*2)
    df.ii = bws.i$SDF@data[,1:4]
  	names(df.ii) = c("b0", "b1", "b2", "b3")
    rsq.i = c(aicval, rsq(gr$b1, df.ii$b1), rsq(gr$b2, df.ii$b2),
              rsq(gr$b3, df.ii$b3))
    rmse.i = c(rmse(gr$b0, df.ii$b0), rmse(gr$b1, df.ii$b1), rmse(gr$b2, df.ii$b2),
               rmse(gr$b3, df.ii$b3))
    mae.i = c(mae(gr$b0, df.ii$b0), mae(gr$b1, df.ii$b1), mae(gr$b2, df.ii$b2),
              mae(gr$b3, df.ii$b3))
  	mgwr_res = c(rsq.i, rmse.i, as.vector(mae.i))
    names(mgwr_res) = c("AIC", "Rsq_B1", "Rsq_B2", "Rsq_B3", "RMSE_B0", "RMSE_B1",
    		"RMSE_B2","RMSE_B3", "MAE_B0", "MAE_B1", "MAE_B2","MAE_B3")
    ## model prediction accuracy
    mgwr.pred = as.vector(bws.i$SDF$yhat)
    rsq.mgwr.y = rsq(mgwr.pred, df.i$y)
    rmse.mgwr.y = rmse(mgwr.pred, df.i$y)
    mae.mgwr.y = mae(mgwr.pred, df.i$y)
    ypred.mgwr = c(rsq.mgwr.y, rmse.mgwr.y, mae.mgwr.y)
    names(ypred.mgwr) = c("R2", "RMSE", "MAE")
    bws = as.vector(bws.i$GW.arguments$bws)
    names(bws) = c("bw_mgwr_b0", "bw_mgwr_b1", "bw_mgwr_b2", "bw_mgwr_b3")

    return(list(mgwr_res = mgwr_res, mgwr_resids = resids, ypred_mgwr = ypred.mgwr, bws = bws))
}

# Loads precompiled results 
# Edit EvalSim = T if you want to create from fresh 
EvalSim = F
if(EvalSim) {
  ## Do GAM
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  gam_list <- foreach (i = 1:100, .packages = c("tidyverse")) %dopar% {
        gam_mods_accuracies(dogam_function(i))
  }
  stopCluster(cl)
  length(gam_list); names(gam_list[[1]])
  save(gam_list, file = "parallel_gam_res_paper2.RData")
  ## Do MGWR
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  mgwr_list <- foreach (i = 1:100) %dopar% {
        mgwr_mods_accuracies(dogwr_function(i))
  }
  stopCluster(cl)
  length(mgwr_list); names(mgwr_list[[1]])
  save(mgwr_list, file = "parallel_mgwr_res_paper2.RData")
  } else {
    load("parallel_gam_res_paper2.RData") # 30 minutes to run
    load("parallel_mgwr_res_paper2.RData") # 10.5 hrs to run
}

# create summary tables
gam_res_tab = matrix(unlist(lapply(gam_list, function(x) (x$gam_res))), nrow = 100, byrow = T)
mgwr_res_tab = matrix(unlist(lapply(mgwr_list, function(x) (x$mgwr_res))), nrow = 100, byrow = T)
colnames(gam_res_tab) = c("AIC", "Rsq_B1", "Rsq_B2", "Rsq_B3", "RMSE_B0", "RMSE_B1", "RMSE_B2","RMSE_B3",
                      "MAE_B0", "MAE_B1", "MAE_B2","MAE_B3")
colnames(mgwr_res_tab) = c("AIC", "Rsq_B1", "Rsq_B2", "Rsq_B3","RMSE_B0", "RMSE_B1", "RMSE_B2", "RMSE_B3",
                           "MAE_B0", "MAE_B1", "MAE_B2", "MAE_B3")
gam_res_tab = data.frame(gam_res_tab, res = "GGP-GAM")
mgwr_res_tab = data.frame(mgwr_res_tab, res = "MGWR")

# resids
gam_resids = matrix(unlist(lapply(gam_list, function(x) (x$gam_resids))), ncol = 100, byrow = F)
mgwr_resids = matrix(unlist(lapply(mgwr_list, function(x) (x$mgwr_resids))), ncol = 100, byrow = F)

# y pred
ypred_gam = matrix(unlist(lapply(gam_list, function(x) (x$ypred_gam))), nrow = 100, byrow = T)
ypred_mgwr = matrix(unlist(lapply(mgwr_list, function(x) (x$ypred_mgwr))), nrow = 100, byrow = T)
# join
ypred = cbind(ypred_gam, ypred_mgwr)
colnames(ypred) = c("R2_gam", "Rmse_gam", "Mae_gam", "R2_mgwr", "Rmse_mgwr", "Mae_mgwr")

# sps bws
sps = matrix(unlist(lapply(gam_list, function(x) (x$sps))), nrow = 100, byrow = T)
bws = matrix(unlist(lapply(mgwr_list, function(x) (x$bws))), nrow = 100, byrow = T)
bws_sp = data.frame(sps, bws)
colnames(bws_sp) = c("gam_b0", "gam_b1", "gam_b2", "gam_b3", "mgw_b0", "mgw_b1", "mgw_b2", "mgw_b3")

## Figure 5
# join results form gam and mgwr
df = rbind(gam_res_tab, mgwr_res_tab)
rownames(df) = 1:nrow(df)

df %>% 
  pivot_longer(-res) -> df_long
## now carry on! 
df_long$name = recode(df_long$name,
  "MAE_B0" = '"MAE "*beta[0]',
  "MAE_B1" = '"MAE "*beta[1]',
  "MAE_B2" = '"MAE "*beta[2]',
  "MAE_B3" = '"MAE "*beta[3]',
  "RMSE_B0" = '"RMSE "*beta[0]',
  "RMSE_B1" = '"RMSE "*beta[1]',
  "RMSE_B2" = '"RMSE "*beta[2]',
  "RMSE_B3" = '"RMSE "*beta[3]',
  "Rsq_B1"= 'R^2*" "*beta[1]',
  "Rsq_B2"= 'R^2*" "*beta[2]',
  "Rsq_B3"= 'R^2*" "*beta[3]')

p4 = 
  df_long %>%   
  #filter(name != "AIC") %>%
  ggplot(aes(y = value,  fill = res)) +
  geom_boxplot(lwd = .2, outliers = FALSE) +
  scale_fill_discrete_c4a_cat(palette="tableau.classic_traffic_light", name = "") + 
  facet_wrap(~name, scales = "free", labeller=label_parsed, dir = "h", as.table = T) + 
  theme_bw() +
  ylab("") +
  theme(#legend.position = c(0.85, 0.15),
       legend.position = "bottom", 
       legend.key.size = unit(1.4, "cm"),
       legend.text=element_text(size=10), 
       axis.title.x = element_blank(),
	     axis.text.x = element_blank(),
       axis.ticks.x = element_blank())
p4

# Table 1
# Investigations of the large differences in AIC values
tab1 = round(rbind(summary(as.vector(gam_resids)), summary(as.vector(mgwr_resids))), 3)
tab1 = format(tab1, digits = 5)
rownames(tab1) = c("GGP-GAM", "MGWR")
knitr::kable(tab1, digits = 3, caption = "\\label{tab:tab1}Summaries of the residuals of the 100 GGP-GAM and MGWR models.", linesep = "") 

# Table 2
# The spatial variation in the residuals was also compared
library(spdep) 
# create point layer
gr_sf = st_as_sf(gr,coords=c("px","py"))
gr_sf$ID = 1: nrow(gr_sf)
# create voronoi polys
gr_pols = st_collection_extract(st_voronoi(do.call(c, st_geometry(gr_sf))))
 # match them to points:
gr_sf$pols = gr_pols[unlist(st_intersects(gr_sf, gr_pols))]
# assign 
gr_pols = st_set_geometry(gr_sf, "pols")["ID"]
# do moran test
gr.nb <- poly2nb(gr_pols)
gr.lw = nb2listw(gr.nb)
gam_moran = as.vector(apply(gam_resids, 2, 
                             function(x) moran.test(x, gr.lw)$estimate))
mgwr_moran = as.vector(apply(mgwr_resids, 2, 
                             function(x) moran.test(x, gr.lw)$estimate))
tab2 = rbind(summary(gam_moran),
            summary(mgwr_moran))
rownames(tab2) = c("GGP-GAM", "MGWR")
negative SA (i.e. perfect dispersion
knitr::kable(tab2, digits = 3, caption = "\\label{tab:tab2}Summaries of the Moran's I of the residuals of the 100 GGP-GAM and MGWR models.", linesep = "") 

## Figure 6
# The accuracy of the predicted measures of $y$ (i.e. $\hat{y}$) arising from the two models 
# seperate for facet plot  
ypred %>%
  data.frame() %>%
  select(ends_with("_gam")) %>% 
  transmute(R2 = R2_gam, RMSE = Rmse_gam, MAE = Mae_gam) -> gamy 
ypred %>%
  data.frame() %>%
  select(ends_with("_mgwr")) %>% 
  transmute(R2 = R2_mgwr, RMSE = Rmse_mgwr, MAE = Mae_mgwr) -> mgwy 
# add label bind and tidy
gamy$res = "GGP-GAM"
mgwy$res = "MGWR"
df = rbind(gamy, mgwy)
rownames(df) = 1:nrow(df)
# make long and recode names
df %>% 
  pivot_longer(-res) -> df_long
df_long$name = recode(df_long$name,
  "MAE" = "MAE",
  "RMSE" = "RMSE",
  "R2"= 'R^2')
# plot
p5 = 
  df_long %>%   
  ggplot(aes(y = value,  fill = res)) +
  geom_boxplot(lwd = .2) +
  scale_fill_discrete_c4a_cat(palette="tableau.classic_traffic_light", name = "") + 
  facet_wrap(~factor(name, levels = c('R^2', "MAE", "RMSE")), 
             scales = "free", labeller=label_parsed, dir = "h", as.table = T) + 
  theme_bw() +
  ylab("") +
  theme(#legend.position = c(0.85, 0.15),
       legend.position = "bottom", 
       legend.key.size = unit(1.4, "cm"),
       legend.text=element_text(size=10), 
       axis.title.x = element_blank(),
	   axis.text.x = element_blank(),
       axis.ticks.x = element_blank())
p5

## Table 3
bws_sp %>% 
  mutate(ID = 1:n()) %>%
  pivot_longer(-ID) -> df
df$name = recode(df$name,
  "gam_b0" = '"GAM SP "*x[0]',
  "gam_b1" = '"GAM SP "*x[1]',
  "gam_b2" = '"GAM SP "*x[2]',
  "gam_b3" = '"GAM SP "*x[3]',
  "mgw_b0" = '"MGWR BW "*x[0]',
  "mgw_b1" = '"MGWR BW "*x[1]',
  "mgw_b2" = '"MGWR BW "*x[2]',
  "mgw_b3" = '"MGWR BW "*x[3]')

# extract the GAM SP and MGWR BW summaries
tab3 = t(apply((bws_sp), 2, summary))
tab3 = data.frame(tab3)
rownames(tab3) = c("GGP-GAM SP $x$~0~", "GGP-GAM SP $x$~1~", "GGP-GAM SP $x$~2~", "GGP-GAM SP $x$~3~", 
                  "MGWR BW $x$~0~", "MGWR BW $x$~1~", "MGWR BW $x$~2~", "MGWR BW $x$~3~")
# format the GAM SPs to scientific notation
tab3[1:4,] = format(tab3[1:4,], digits = 2, scientific = T, big.mark = ",")
# for the mean
tab3 = tab3[,-4]
colnames(tab3) = c("Min.", "Q1", "Median", "Q3", "Max.")
knitr::kable(tab3, digits = 3, caption = "\\label{tab:tab3}Summaries of MGWR bandwidths (BW) and GGP-GAM spline smoothing parameters (SP).", linesep = "") 

#' ## 4.2 A single GGP-GAM in detail
#' 
## Table 4
load("sim_paper2.RData")
i = 51
df.i = data.frame(y = y_sim[,i],
                      x1 = x1_sim[,i],
                      x2 = x2_sim[,i],
                      x3 = x3_sim[,i],
                      Intercept = 1, 
                      X = gr$px, Y = gr$py)
kval = 155
gam51 <- gam(y~ 0 + 
                Intercept + s(X,Y,bs='gp',k=kval,by=Intercept) + 
                x1 + s(X,Y,bs='gp',k=kval,by=x1) +
                x2 + s(X,Y,bs='gp',k=kval,by=x2) + 
                x3 + s(X,Y,bs='gp',k=kval,by=x3), 
                data = df.i, 
                method = "REML")
terms = c("Intercept", "x1", "x2", "x3")
df51_svc = ggp_gam_svc(gam51, terms, df.i)
tab4 = as.data.frame(k.check(gam51, subsample = 5000, n.rep = 200))
knitr::kable(tab4, digits = 3, caption = "\\label{tab:tab4}Diagnostics of the GGP-GAM smoothing optimisation for the SVCs of a single GGP-GAM model.", linesep = "", booktabs = T) 

## Table 5
tab5 = as.data.frame(round(summary(gam51)$p.table, 3))
knitr::kable(tab5, digits = 3, caption = "\\label{tab:tab5}The fixed (parametric) coefficients for a single GGP-GAM.", linesep = "", booktabs = T) 

## Table 6
tab6 = as.data.frame(round(summary(gam51)$s.table, 3))
knitr::kable(tab6, digits = 3, caption = "\\label{tab:tab6}The smooth terms of the GGP splines for a single GGP-GAM model.", linesep = "", booktabs = T) 

## Table 7
df51_svc %>% 
  select(starts_with("b_")) %>% 
  sapply(summary) %>% 
  t() -> tab7
rownames(tab7) = c("$\\beta_{0}$", "$\\beta_{x1}$", "$\\beta_{x2}$", "$\\beta_{x3}$")
knitr::kable(tab7, digits = 3, caption = "\\label{tab:tab7}Summaries of the spatially varying coefficients for a single GGP-GAM model.", linesep = "", booktabs = T) 

## Figure 7
# join gam svc
gr = cbind(gr, gam_b0 = df51_svc$b_Intercept)
gr = cbind(gr, gam_b1 = df51_svc$b_x1)
gr = cbind(gr, gam_b2 = df51_svc$b_x2)
gr = cbind(gr, gam_b3 = df51_svc$b_x3)

# Create mgwr svc
# Loads precompiled results 
# Edit DoAnal = T if you want to create from fresh 
DoAnal = FALSE
if(DoAnal) {
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  mgwr_mod <- foreach (i = 51) %dopar% {
        dogwr_function(i)
    }
  stopCluster(cl)
  save(mgwr_mod, file = "mgwr_mod.RData")
  } else {
    load("mgwr_mod.RData")
} 

# join mgwr svc
# head(mgwr_mod[[1]]$SDF)
gr = cbind(gr, mgw_b0 = mgwr_mod[[1]]$SDF$Intercept)
gr = cbind(gr, mgw_b1 = mgwr_mod[[1]]$SDF$x1)
gr = cbind(gr, mgw_b2 = mgwr_mod[[1]]$SDF$x2)
gr = cbind(gr, mgw_b3 = mgwr_mod[[1]]$SDF$x3)

# create the GGP-GAM plots 
tit =expression(paste("GGP-GAM "*beta[0]*""))
p6 = beta_plot_function(gr, beta = gr$gam_b0, tit) 
tit =expression(paste("GGP-GAM "*beta[1]*""))
p7 = beta_plot_function(gr, beta = gr$gam_b1, tit) 
tit =expression(paste("GGP-GAM "*beta[2]*""))
p8 = beta_plot_function(gr, beta = gr$gam_b2, tit) 
tit =expression(paste("GGP-GAM "*beta[3]*""))
p9 = beta_plot_function(gr, beta = gr$gam_b3, tit) 

# create the MGWR plots 
tit =expression(paste("MGWR "*beta[0]*""))
p10 = beta_plot_function(gr, beta = gr$mgw_b0, tit) 
tit =expression(paste("MGWR "*beta[1]*""))
p11 = beta_plot_function(gr, beta = gr$mgw_b1, tit) 
tit =expression(paste("MGWR "*beta[2]*""))
p12 = beta_plot_function(gr, beta = gr$mgw_b2, tit) 
tit =expression(paste("MGWR "*beta[3]*""))
p13 = beta_plot_function(gr, beta = gr$mgw_b3, tit) 

plot_grid(p0, p1, p2, p3, legend_betas,
          p6, p7, p8, p9, legend_betas, 
          p10, p11, p12, p13, legend_betas, nrow = 3, rel_widths = c(1,1,1,1,0.25))

#' ## 4.3 GGP-GAM tuning with a larger datase

## set up dataset
i = 1
df.i = data.frame(y = sims$y_sim[,i],
                      x1 = sims$x1_sim[,i],
                      x2 = sims$x2_sim[,i],
                      x3 = sims$x3_sim[,i],
                      Intercept = 1, 
                      X = sims$px, Y = sims$py)

# This takes 13 hours to run 
# Loads precompiled results 
# Edit DoBigBam = T if you want to create from fresh 
DoBigBam = F
if(DoBigBam) {
  klist = c(100, 250, 500, 750, 1000, 1500, 2000)
  for (i in 1:length(klist)){
  	kval = klist[i]
  	cl <- makeCluster(detectCores()-1)
  	ti = Sys.time()
  	gam.ii <- bam(y~ 0 + 
                  Intercept + s(X,Y,k = kval,bs='gp',by=Intercept) + 
                  x1 + s(X,Y,k = kval,bs='gp',by=x1) +
                  x2 + s(X,Y,k = kval,bs='gp',by=x2) + 
                  x3 + s(X,Y,k = kval,bs='gp',by=x3), 
                  data = df.i, 
                  method = "fREML",cluster=cl)
  	tfin = Sys.time() - ti 
  	stopCluster(cl)
  	tit = paste0("gam.ii.", kval, ".RData")
  	save(tfin, gam.ii, file = tit)
  	cat(i, "\t")
  }
}   

# Loads precompiled results 
# Edit DoBigModelAnalysis = T if you want to create from fresh 
# BUT you will need to create the .RData files using the code above (DoBigBam)
# This is because these are too big to share 
# from 28MB (gam.ii.100.RData) to 2.13GB (gam.ii.2000.RData)
DoBigModelAnalysis = F
if(DoBigModelAnalysis) {
  # load the GGP-GAMs
  load("gam.ii.100.RData"); gam100 <- gam.ii; t100 <- tfin
  load("gam.ii.250.RData"); gam250 <- gam.ii; t250 <- tfin
  load("gam.ii.500.RData"); gam500 <- gam.ii; t500 <- tfin
  load("gam.ii.750.RData"); gam750 <- gam.ii; t750 <- tfin
  load("gam.ii.1000.RData"); gam1000 <- gam.ii; t1000 <- tfin
  load("gam.ii.1500.RData"); gam1500 <- gam.ii; t1500 <- tfin
  load("gam.ii.2000.RData"); gam2000 <- gam.ii; t2000 <- tfin
  timings = as.vector(c(t100, t250, t500, t750, t1000, t1500, t2000)/60)
  k = c(100, 250, 500, 750, 1000, 1500, 2000)
  gam.list = c("gam100", "gam250", "gam500", "gam750", "gam1000", "gam1500", "gam2000")
  # model measures
  singles = matrix(0,nrow = length(gam.list),ncol = 5)
  colnames(singles) = c("deviance", "null.deviance", "df.residual", "rank", "sig2")
  rownames(singles) = k
  
  resi = fits = sps = vector()
  edfs = param_coef.est = param_coef.pval = sp_pvals = vector()
  
  for(i in 1:length(gam.list)) {
    gam.i = get(gam.list[i])
    fitted.vals = gam.i$fitted.values
    resi = cbind(resi, gam.i$residuals)
    # model characteristics
    singles[i,] = 
      c(gam.i$deviance, #model deviance (not penalized deviance).
    gam.i$null.deviance, 	#deviance for single parameter model.
    gam.i$df.residual, #effective residual degrees of freedom of the model.
    gam.i$rank,  #apparent rank of fitted model.
    gam.i$sig2)	# estimated or supplied variance/scale parameter.
    # converged = gam.i$converged T/F
    sps = rbind(sps, gam.i$sp)  # smoothing parameters
    # fits
    resids = gam.i$residuals
    aicval = nrow(df.i) * (log(2*pi)+1 + log( (sum(resids^2)/nrow(df.i) ) ) ) + ((4+1)*2)
    r2val = rsq(df.i$y, fitted.vals)
    rmseval = rmse(df.i$y, fitted.vals)
    maeval = mae(df.i$y, fitted.vals)
    fits = rbind(fits, c(aicval, r2val,rmseval,maeval))
    
    # model checks
    edfs = rbind(edfs, as.data.frame(k.check(gam.i, subsample = 5000, n.rep = 200))$edf)
    param_coef.est = rbind(param_coef.est, as.data.frame(summary(gam.i)$p.table)$Estimate)
    param_coef.pval = rbind(param_coef.pval, as.data.frame(summary(gam.i)$p.table)$`Pr(>|t|)`)
    sp_pvals = rbind(sp_pvals, as.data.frame(summary(gam.i)$s.table)$`p-value`)    
  }
  # do naming
  colnames(resi) = gam.list
  colnames(sps) = paste0("splineSP_", c("x0", "x1", "x2", "x3") )
  colnames(fits) = c("AIC", "R2", "RMSE", "MAE")
  rownames(sps) = rownames(fits) = k
  
  rownames(edfs) = rownames(param_coef.est) = rownames(param_coef.pval) = rownames(sp_pvals) = k
  colnames(param_coef.est)  = paste0("paramCoef_", c("x0", "x1", "x2", "x3")) 
  colnames(param_coef.pval) = paste0("paramPval_", c("x0", "x1", "x2", "x3")) 
  colnames(edfs) = paste0(c("spXYx0", "spXYx1", "spXYx2", "spXYx3"), "_edf")
  colnames(sp_pvals) = paste0(c("spXYx0", "spXYx1", "spXYx2", "spXYx3"), "_pval")
  
  # bind together and save
  df = data.frame(k = k, time_mins = timings, singles, fits, sps,
                  param_coef.est, param_coef.pval, edfs, sp_pvals, row.names = NULL)
  save(df, resi, file = "bigmodelanalysis.RData")
} else {
  load("bigmodelanalysis.RData")
}

## Figure 8
equal_breaks <- function(n = 3, s = 0.05, r = 0,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    if(seq[2]-seq[1] < 10^(-r)) seq else round(seq, r)
  }
}
formatter <- function(...){
  function(x) format(round(x, 1), ...)
}
df %>% transmute(k, 
                 "Time (mins)"= as.vector(time_mins), 
                 "Model Deviance" = deviance, 
                 "Residual DoF" = df.residual, 
                 "Estimated Variance" = round(sig2,5), 
                 AIC, 
                 R2 = round(R2,5), 
                 RMSE = round(RMSE,5), 
                 MAE = round(MAE,5)) %>% 
  pivot_longer(-k) %>% 
  ggplot(aes(x = value, y = k)) +
  geom_point() +
  geom_line() +
  facet_wrap(~factor(name, levels = c("Time (mins)", "Model Deviance", 
                                      "Residual DoF", "Estimated Variance", 
                                      "AIC", "R2", "RMSE", "MAE")), 
             scales = "free", ncol = 4) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
  xlab("") + 
  theme_bw()

## Table 8
tab8 = t(apply(resi, 2, function(x) round(summary(x),4)))
tab8 = cbind(k = c(100, 250, 500, 750, 1000, 1500, 2000), tab8)
knitr::kable(tab8,caption = "\\label{tab:tab8}Summaries of the GGP-GAM residuals with increasing values of k.",
             linesep = "", booktabs = T,  row.names = F) 

## Table 9
df %>% select(k, starts_with("param")) %>% 
  relocate(k,
           x0 = paramCoef_x0, `x0 p-val` = paramPval_x0, 
           x1 = paramCoef_x1, `x1 p-val` = paramPval_x1,
           x2 = paramCoef_x2, `x2 p-val` = paramPval_x2, 
           x3 = paramCoef_x3, `x3 p-val` = paramPval_x3) %>% 
  round(3) |> 
  
  mutate(k = as.character(k)) %>%
  mutate_if(is.numeric, format, digits=3,nsmall = 3) %>%
  knitr::kable(digits = 3, caption = "\\label{tab:tab9}The GGP-GAM parametric coefficients and their global significance (p-val) with increasing values of k.", linesep = "", booktabs = T, row.names = F) 

## Table 10
df %>% select(k, starts_with("spXY")) %>%  
  relocate(k,
           x0 = spXYx0_edf, `x0 p-val` = spXYx0_pval, 
           x1 = spXYx1_edf, `x1 p-val` = spXYx1_pval,
           x2 = spXYx2_edf, `x2 p-val` = spXYx2_pval, 
           x3 = spXYx3_edf, `x3 p-val` = spXYx3_pval) %>%
  round(3) %>%
  mutate(k = as.character(k)) %>%
  mutate_if(is.numeric, format, digits=3,nsmall = 3) %>%
  knitr::kable(digits = 3, caption = "\\label{tab:tab10}The effective degrees of freedom of the GGP spline smooth terms and their local significance (p-val) with increasing values of k.", linesep = "", booktabs = T, row.names = F) 

## Table 11
df %>% select(k, starts_with("spline")) -> tab11
names(tab11) = c("k", rownames(tab5))
tab11 %>% 
  #round(3) %>%
  mutate(k = as.character(k)) %>%
  mutate_if(is.numeric, format, digits=3,nsmall = 3) %>%
  knitr::kable(digits = 3, caption = "\\label{tab:tab11}The GGP spline smoothing parameters with increasing k.", linesep = "", booktabs = T, row.names = F) 

#' ## 4.4 Empirical example: Brexit vote

library(stgam)
library(GWmodel)
brex_gam = gam(share_leave ~ 0 + one +
                 s(X, Y, bs='gp',by=one) +
                 s(X, Y, bs='gp',by=christian) +
                 s(X, Y, bs='gp',by=degree_educated) +
                 s(X, Y, bs='gp',by=no_car) +
                 s(X, Y, bs='gp',by=younger_adults),
               data = gb_df, method = "REML")
# calculate the SVCs using stgam function
terms = c("one","christian", "degree_educated", "no_car", "younger_adults")
svcs = calculate_vcs(gb_df, brex_gam, terms)

# MGWR - this will take a few minutes to run
# convert to sp
gb.sp <- gb |> as("Spatial")
# calculate bandwidths
bws.i = gwr.multiscale(share_leave~christian+degree_educated+no_car+younger_adults,  
  	                   data = gb.sp,
                       kernel = "bisquare",
  	                   bws0=c(1000,1000,1000,1000,1000))
msgwr_bws.i = as.vector(bws.i$GW.arguments$bws)
# apply bandwidths to model
brex_mgwr = gwr.multiscale(share_leave~christian+degree_educated+no_car+younger_adults,
  	                       data = gb.sp,
  	                       kernel = "bisquare",
                           bws0 = msgwr_bws.i,
                           bw.seled=rep(T, 5))

## Table 12
r2 = c(rsq(brex_mgwr$SDF$yhat, gb.sp$share_leave), rsq(brex_gam$fitted.values,gb_df$share_leave))
m = c(mae(brex_mgwr$SDF$yhat, gb.sp$share_leave), mae(brex_gam$fitted.values,gb_df$share_leave))
r = c(rmse(brex_mgwr$SDF$yhat, gb.sp$share_leave), rmse(brex_gam$fitted.values,gb_df$share_leave))
a = c(nrow(gb.sp) * (log(2*pi)+1 + log( (sum(brex_mgwr$SDF$residual^2)/nrow(gb.sp) ) ) ) + ((4+1)*2), 
      nrow(gb_df) * (log(2*pi)+1 + log( (sum(brex_gam$residuals^2)/nrow(gb_df) ) ) ) + ((4+1)*2) )
tab12 = data.frame(R2 = round(r2, 3), MAE = round(m, 3), RMSE = round(r, 3), AIC = round(a, 1), 
                row.names = c("MGWR", "GGP-GAM"))    
tab12 |>
  knitr::kable(caption = "\\label{tab:tab12}Fit measures from the GGP-GAM and MGWR models of Brexit vote.", linesep = "", booktabs = T, row.names = T) 

# create GAM SVCs
library(stgam) 
terms = c("one","christian", "degree_educated", "no_car", "younger_adults")
gam_svcs = calculate_vcs(gb_df, brex_gam, terms)

gb_gam <- 
  gb |> select(lad15cd) |> 
  left_join(gam_svcs |> select(lad15cd, all_of(starts_with('b_')))) 
gb_mgwr <- 
  gb |> select(lad15cd) |> 
  left_join(data.frame(lad15cd = gb.sp$lad15cd, brex_mgwr$SDF@data[, 1:5]))

# local GAM p-values and MGWR bandwidths
gam_pvals = summary(brex_gam)[[8]]
mgwr_bws = msgwr_bws.i/1000

## Table 13
gb_gam %>% st_drop_geometry() |>
  select(starts_with("b_")) |> 
  sapply(summary) |> 
  rbind(`Smooth p-values` = gam_pvals) |>
  #rbind(`SVC range (km)` = round(gam_dists, 1)) |>
  t() -> tab13
rownames(tab13) = c("$\\beta_{0}$", "$\\beta_{Christian}$", "$\\beta_{Degree}$", "$\\beta_{No Car}$", "$\\beta_{Younger}$")
knitr::kable(tab13, digits = 3, caption = "\\label{tab:tab13}Summaries of the SVCs for the Brexit GGP-GAM model.", linesep = "", booktabs = T) 

## Table 14
gb_mgwr %>% st_drop_geometry() |>
  select(Intercept:younger_adults) |>
  sapply(summary) |>
  rbind(`Bandwidth (km)` = round(mgwr_bws,1)) |>
  #rbind(`SVC range (km)` = round(mgwr_dists, 1)) |>
  t() -> tab14
rownames(tab14) = c("$\\beta_{0}$", "$\\beta_{Christian}$", "$\\beta_{Degree}$", "$\\beta_{No Car}$", "$\\beta_{Younger}$")
knitr::kable(tab14, digits = 3, caption = "\\label{tab:tab14}Summaries of the SVCs for the Brexit MGWR model.", linesep = "", booktabs = T) 

## Figure 9
svc_plot_function = function(gb_sf = gb_gam, beta = gb_gam$b_one, tit, legend.pos = "bottom",
                             c4apal = "kovesi.bu_wh_rd2") {
  ggplot(gb_sf, aes(fill = beta)) + geom_sf(col = NA) +
  scale_fill_continuous_c4a_div(palette=c4apal, mid = 0) + 
  theme_bw() + xlab("") + ylab("") +  
  labs(subtitle = tit) +
  theme(legend.position = "none", 
        legend.title = element_blank(),
        legend.key.width = unit(.5, "cm"),
        legend.text = element_blank(),
        legend.key.height = unit(.35, "cm"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 
}

tit = expression(paste("GGP-GAM "*beta[0]*""))
pg_0 = svc_plot_function(gb_gam, beta = gb_gam$b_one, tit) 
tit = expression(paste("GGP-GAM "*beta[Christian]*""))
pg_1 = svc_plot_function(gb_gam, beta = gb_gam$b_christian, tit) 
tit = expression(paste("GGP-GAM "*beta[Degree]*""))
pg_2 = svc_plot_function(gb_gam, beta = gb_gam$b_degree_educated, tit) 
tit = expression(paste("GGP-GAM "*beta[`No Car`]*""))
pg_3 = svc_plot_function(gb_gam, beta = gb_gam$b_no_car, tit) 
tit = expression(paste("GGP-GAM "*beta[`Younger`]*""))
pg_4 = svc_plot_function(gb_gam, beta = gb_gam$b_younger_adults, tit) 

tit = expression(paste("MGWR "*beta[0]*""))
pm_0 = svc_plot_function(gb_mgwr, beta = gb_mgwr$Intercept, tit) 
tit = expression(paste("MGWR "*beta[Christian]*""))
pm_1 = svc_plot_function(gb_mgwr, beta = gb_mgwr$christian, tit) 
tit = expression(paste("MGWR "*beta[Degree]*""))
pm_2 = svc_plot_function(gb_mgwr, beta = gb_mgwr$degree_educated, tit) 
tit = expression(paste("MGWR "*beta[`No Car`]*""))
pm_3 = svc_plot_function(gb_mgwr, beta = gb_mgwr$no_car, tit) 
tit = expression(paste("MGWR "*beta[`Younger`]*""))
pm_4 = svc_plot_function(gb_mgwr, beta = gb_mgwr$younger_adults, tit) 

plot_grid(pg_0, pg_1, pg_2, pg_3, pg_4, 
          pm_0, pm_1, pm_2, pm_3, pm_4, 
          ncol = 5)


### END ###