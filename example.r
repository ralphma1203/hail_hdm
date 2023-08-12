source("hail_function.r")

load("data.RData") # load dataset


# create grouping variable (event) from raw data frame 
data_frame <- dat %>% ungroup() %>% mutate(row_ind = row_number()) %>%
	group_by(evt_date, N_COUNTY) %>%
	mutate(event_ind = cur_group_id())


load("pairs_1.RData") # load pairs index

pair.dat <- as.data.frame(pair.dat) # convert to data frame layout

names(pair.dat) <- c( "ind1", "ind2", "distance")

logic.vector <- data.frame %>% ungroup() %>% select(evt_date) %>%
		mutate(year = lubridate::year(evt_date)) %>% select(year)

cat(NROW(data.frame), NROW(pair.dat),"\n")

# Drop 2015 and use list assign here
c(
data.frame,
pair.dat,
event_contribute)	%<-% 
drop.by.vector(data.frame, pair.dat, logic.vector, drop.label=2015)  # drop with evt_date in 2015

cat(NROW(data.frame), NROW(pair.dat),"\n")

form <- as.formula(
	claim ~ hail_size_pif + log10(deductible) + log10(BLDG_LMT) + 
			I(NUM_FAMILIES > 0) + I(YR_BUILT >= 1990) + I(TYP_CONST_CD == "1 Frame") +
			I(ROOF_TYPE_CD == "A Asphalt shingles")
	)
	
unique.event <- unique(data.frame$event_ind)

# Example of dropping one event
c(dat_drop_one_event,
  pair.dat_drop_one_event,
  event_contribute) %<-% drop.one.event(data.frame, pair.dat, event_index=unique.event[ind])

cat(NROW(dat_drop_one_event), NROW(pair.dat_drop_one_event),"\n")	

rm(data.frame, pair.dat)

load("t14_12.RData")

ini <- tail(ans$est,3)

t0 <- proc.time()

# refitting model and storing gradient and hessian matrix (optional)
ans <- est1.twostage(dat_drop_one_event, pair.dat_drop_one_event, form, grad=TRUE, hessian=TRUE,
		ini, tol.bound=1e-5, factr=1e8, progress=TRUE)
			
compute.time <- proc.time() - t0

radius <- 1

save(ans, compute.time, event_contribute, form, radius, file=paste0("t14_12_e_", unique.event[ind],".RData"))


#################### Prediction example ####################
load("data.RData") # load dataset

dat.claim <- dat4 %>% ungroup() %>%  mutate(total_loss = loss_tot_amnt_paid + deductible) %>% 
  filter(claim == 1 & loss_tot_amnt_paid > 0) %>% 
  mutate(row_ind = row_number()) %>%
  group_by(evt_date, N_COUNTY) %>%
  mutate(event_ind = cur_group_id())

load("CO_claim1_loss_positive_pairs_100.RData") # load pairs

pair.dat.full <- as.data.frame(pair.dat) # convert to data frame layout

names(pair.dat.full) <- c( "ind1", "ind2", "distance")

ind.event <- 79 # event ID 79: 1188 records
radius <- 50

ind.event <- 115 # 36 records

cat(ind.event, radius, "\n")



pair.dat <- pair.dat.full %>% filter(distance <= radius)


logic.vector <- dat.claim %>% ungroup() %>% select(evt_date) %>%
		mutate(year = lubridate::year(evt_date)) %>% select(year)

cat(NROW(dat.claim), NROW(pair.dat),"\n")

# Drop 2015 and use list assign here
c(
dat.claim,
pair.dat,
event_contribute)	%<-% 
drop.by.vector(dat.claim, pair.dat, logic.vector, drop.label=2015)  # drop with evt_date in 2015

cat(NROW(dat.claim), NROW(pair.dat),"\n")


form <- as.formula(
	Surv(days_to_rpt2, rep(1,NROW(dat_drop_one_event))) ~ 
			hail_size_pif + log10(deductible) + log10(BLDG_LMT) + # "obvious" relationship
			I(NUM_FAMILIES > 0) + I(YR_BUILT >= 1990) + I(TYP_CONST_CD == "1 Frame") +
			I(ROOF_TYPE_CD == "A Asphalt shingles") # seem not useful
	)

unique.event <- unique(dat.claim$event_ind)
dat.claim %>% group_by(event_ind) %>% tally() %>% arrange(desc(n)) %>% ungroup() %>% head(50)

# Just pick one event for checking
ind.event <- 79 # event ID 79: 1188 records
ind.event <- 115 # 36 records

c(dat_drop_one_event,
  pair.dat_drop_one_event,
  event_contribute) %<-% drop.one.event(dat.claim, 
	pair.dat, event_index=ind.event)
	
c(dat_single_event,
  pair.dat_single_event) %<-% single.event(dat.claim, 
	pair.dat, event_index=ind.event)	

cat(NROW(dat_drop_one_event), NROW(pair.dat_drop_one_event),"\n")
cat(NROW(dat_single_event), NROW(pair.dat_single_event),"\n")
  

# prepare design matrix and response
form.single <- as.formula(
	Surv(days_to_rpt2, rep(1,NROW(dat_single_event))) ~ 
			hail_size_pif + log10(deductible) + log10(BLDG_LMT) + # "obvious" relationship
			I(NUM_FAMILIES > 0) + I(YR_BUILT >= 1990) + I(TYP_CONST_CD == "1 Frame") +
			I(ROOF_TYPE_CD == "A Asphalt shingles") # seem not useful
	)


obj_tmp <- survival::survreg(form.single, data=dat_single_event, 
	dist="weibull",x=TRUE, y=TRUE)

	
X_event <- as.matrix(model.matrix(obj_tmp))
y_event <- as.matrix(obj_tmp$y)[,1]
evt_ind <- unique(dat_single_event$event_ind)
within_dist <- 10



longlat <- dat_single_event %>% ungroup() %>% select(LON,LAT)

r_max <- 5
n_node <- 50 ; alpha <- 2
fun=geosphere::distGeo;  progress=TRUE


para_est <- ans_full[[radius]]$est # suppose we have the parameter estimate given radius parameter


df_tmp <- dat_single_event %>% 
	select(days_to_rpt2, loss_tot_amnt_paid,LON, LAT)

pred_2 <- pred2.nnwithin.single.event(X_event, evt_ind, y_event, within_dist=1,  
		para_est, longlat, r_max=8,
		n_node=50, alpha=2,
		fun=geosphere::distGeo,  pred_density=TRUE, pred_seq=seq(1e-2,1e3, by=1e-2), progress=TRUE)
		
##### k = 3  #####
para_est <- ans_full[[radius]]$est


form_T <- as.formula(
	Surv(days_to_rpt2, rep(1,NROW(dat_single_event))) ~ 
			hail_size_pif + log10(deductible) + log10(BLDG_LMT) + # "obvious" relationship
			I(NUM_FAMILIES > 0) + I(YR_BUILT >= 1990) + I(TYP_CONST_CD == "1 Frame") +
			I(ROOF_TYPE_CD == "A Asphalt shingles") # seem not useful
	)

form_Z <- as.formula(
	total_loss ~ 
			hail_size_pif + log10(deductible) + log10(BLDG_LMT) + # "obvious" relationship
			I(NUM_FAMILIES > 0) + I(YR_BUILT >= 1990) + I(TYP_CONST_CD == "1 Frame") +
			I(ROOF_TYPE_CD == "A Asphalt shingles") # seem not useful
	)
	


obj_tmp <- survival::survreg(form_T, data=dat_single_event, 
	dist="weibull",x=TRUE, y=TRUE)

	
X_event_T <- as.matrix(model.matrix(obj_tmp))
T_event <- as.matrix(obj_tmp$y)[,1]


obj_tmp <- glm(form_Z, dat=dat_single_event, family = Gamma(link="log"))

X_event_Z <- as.matrix(model.matrix(obj_tmp))
Z_event <-  obj_tmp$y


evt_ind <- unique(dat_single_event$event_ind)
within_dist <- 10		



pred_3 <- pred3.nnwithin.single.event(
		X_event_T, X_event_Z, T_event, Z_event,
		evt_ind=evt_ind, within_dist=within_dist,  
		para_est=para_est, longlat=longlat, r_max=5,
		in_lower=0, in_upper=1e7, no_subdivision=1e3,
		fun=geosphere::distGeo, pred_density=TRUE, pred_seq=seq(1,1e5, by=1), progress=TRUE)
		
		
		