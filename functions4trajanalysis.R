#==============================================================================
# function to grab all Vineland, ADOS, and Mullen data from specific group
grabGroupData <- function(D,submask,colLabels,grpLabel) {
	# D = data frame with all data read in from main text file
	# submask = mask of subjects to pull out of D
	# colLabels = column labels of interesting data to extract
	# grpLabel = a group label to assign to the extracted data

	varNames = colnames(D)

	tmp_data = D[submask,]
	# correct for the fact that some people have Dx_Count > 5
	tmp_data$Dx_Count[tmp_data$Dx_Count>5] = 5

	# make empty data frames
	data = data.frame(matrix(nrow = sum(tmp_data$Dx_Count), ncol = length(colLabels)))
	colnames(data) = colLabels
	subinfo = data.frame(matrix(nrow = sum(tmp_data$Dx_Count), ncol = 3))
	colnames(subinfo) = c("subjectId","p2f2","Dx")

	# initialize row counter
	row_count = 0
	# loop over subjects
	for (isub in 1:dim(tmp_data)[1]) {
		# get the number of timepoints for current subject
		nTimePoints = tmp_data$Dx_Count[isub]
		# loop over timepoints
		for (itp in 1:nTimePoints) {
			# iterate row_count
			row_count = row_count + 1
			# construct column labels to extract
			labels2use = vector("character", length = 0)
			for (ilabel in 1:length(colLabels)) {
				labels2use = c(labels2use, sprintf("%s_%d",colLabels[ilabel],itp))
			}# for (ilabel in 1:length(allLabels))
			# get a mask of the columns to extract
			colmask = is.element(varNames,labels2use)
			# extract data and store in data frames
			data[row_count,] = tmp_data[isub,colmask]
			rownames(data[row_count,]) = tmp_data$p2f2[isub]
			subinfo[row_count,1:2] = tmp_data[isub,1:2]
			subinfo[row_count,3] = grpLabel
		}# for (itp in 1:nTimePoints)
	}# for (isub in 1:dim(asd_tmp)[1])
	# the final data frame to output
	result = cbind(subinfo,data)
}# function grabGroupData
#==============================================================================



#==============================================================================
# function to clean up any age errors (e.g., age < or = 0)
cleanAgeErrors <- function(df) {
	# df = data frame to process

	# convert age <0 into NA for Vineland
	vals2clean = df$VageMo<=0 & !is.na(df$VageMo)
	df$VageMo[vals2clean] = NA
	# convert age <0 into NA for ADOS
	vals2clean = df$AageMo<=0 & !is.na(df$AageMo)
	df$AageMo[vals2clean] = NA
	# convert age <0 into NA for Mullen
	vals2clean = df$MageMo<=0 & !is.na(df$MageMo)
	df$MageMo[vals2clean] = NA
	df
}# function cleanAgeErrors
#==============================================================================



#==============================================================================
# function to compute lme and make spaghetti plot
spaghettiPlot <- function(df, x_var, y_var, subgrp_var, xLabel, yLabel,
	modelType = "quadratic",
	fname2save = NULL, plot_dots = TRUE, plot_lines = TRUE,
	ci_band = TRUE, pi_band = FALSE, dot_alpha = 1/10, line_alpha = 1/10,
	band_alpha = 1/2, xLimits = c(-5, 60), yLimits = NULL, standardize = FALSE) {
	# DESCRIPTION
	# This function takes a data frame, and strings denoting the x, y, and
	# subgroup variables to plot.  It also computes the linear mixed effect
	# model on the data.
	#
	# INPUT ARGUMENTS
	# df = data frame to process
	# x_var = string denoting variable to plot on x-axis
	# y_var = string denoting variable to plot on y-axis
	# subgrp_var = string denoting variable that has subgroup labels
	# xLabel = label to put on x-axis
	# yLabel = label to put on y-axis
	# fname2save = file path and filename if you want to save plot
	# plot_dots = set to TRUE to plot individual dots
	# plot_lines = set to TRUE to plot individual lines
	# ci_band = set to true if you want to plot confidence bands
	# pi_band = set to true if you want to plot prediction bands
	# dot_alpha = alpha setting for plotting dots
	# line_alpha = alpha setting for plotting lines
	# band_alpha = alpha setting for plotting bands
	# xLimits = set limits on x-axis
	# yLimits = set limits on y-axis
	# standardize = set to TRUE if you want to z-score data before modeling
	#
	# OUTPUT
	# a list with the plot information p and the linear mixed effect model
	#

	# make sure the required libraries are loaded
	require(ggplot2)
	require(nlme)

	# if you are plotting prediction bands, turn off confidence bands
	if (pi_band) {
		ci_band = FALSE
	}# if (pi_band)

	# find unique subgroups and grab them from input data frame
	subgrps = df[,subgrp_var]
	unique_subgrps = unique(subgrps)
	unique_subgrps = unique_subgrps[!unique_subgrps=="NA"]
	#unique_subgrps = unique_subgrps$subgrp[!unique_subgrps$subgrp=="NA"]
	subgrp_mask = is.element(df[,subgrp_var],unique_subgrps)
	df2use = df[subgrp_mask,]
	#df2use = subset(df,is.element(df[[subgrp_var]],unique_subgrps))
	# uniqueSubs = unique(df2use$subjectId)

	#------------------------------------------------------------------------------
	# initialize the plot
	p = ggplot(data = df2use, aes(x = get(x_var), y = get(y_var), group = subjectId))

	#------------------------------------------------------------------------------
	# plot individual lines
	if (plot_lines) {
		p = p + geom_line(aes(colour = get(subgrp_var)), alpha = line_alpha) + guides(alpha = FALSE)
	}# if (plot_lines)

	#------------------------------------------------------------------------------
	# plot individual dots
	if (plot_dots) {
		# plot each data point as a dot
		# p = p + geom_point(aes(colour = get(subgrp_var), alpha = dot_alpha)) + guides(alpha=FALSE)
		p = p + geom_point(aes(colour = get(subgrp_var)), alpha = dot_alpha) + guides(alpha=FALSE)
	} # if (plot_dots)

	#------------------------------------------------------------------------------

	if (modelType=="linear") {
		# compute linear mixed effect model--------------------------------------------
		fx_form = as.formula(sprintf("%s ~ %s*%s", y_var, x_var, subgrp_var))
		rx_form = as.formula(sprintf("~ %s|subjectId",x_var))
	} else if (modelType=="quadratic") {
		# compute quadratic mixed effect model----------------------------------------
		fx_form = as.formula(sprintf("%s ~ %s + I(%s^2) + %s + %s:%s + I(%s^2):%s",
			y_var, x_var, x_var, subgrp_var, x_var, subgrp_var, x_var, subgrp_var))
		# rx_form = as.formula(sprintf("~ %s|subjectId + I(%s^2)|subjectId",x_var, x_var))
		# rx_form = as.formula(sprintf("~ I(%s^2)|subjectId",x_var))
		rx_form = as.formula(sprintf("~ %s|subjectId",x_var))
	}# if (modelType=="linear")

	ctrl <- lmeControl(opt='optim', msMaxIter = 500)
	ml_method = "ML"
	m_allgrps <- eval(substitute(lme(fixed = fx_form, random = rx_form, data = df2use,
		na.action = na.omit, control=ctrl, method = ml_method), list(fx_form = fx_form, rx_form = rx_form)))
	# summary(m_allgrps)
	# anova(m_allgrps)

	#------------------------------------------------------------------------------
	# get information for confidence and prediction bands
	newdat = expand.grid(x = min(df2use[,x_var],na.rm = TRUE):max(df2use[,x_var],na.rm = TRUE), s=sort(unique_subgrps))
	colnames(newdat)[1] = x_var
	colnames(newdat)[2] = subgrp_var
	newdat[,2] = as.factor(newdat[,2])
	newdat$pred <- predict(m_allgrps, newdat, level = 0)
	colnames(newdat)[3] = y_var
	Designmat <- model.matrix(formula(m_allgrps)[-2], newdat)
	predvar <- diag(Designmat %*% vcov(m_allgrps) %*% t(Designmat))
	newdat$SE <- sqrt(predvar)
	newdat$SE2 <- sqrt(predvar+m_allgrps$sigma^2)
	# newdat$subjectId = newdat$subgrpDx
	newdat$subjectId = newdat[,2]

	#------------------------------------------------------------------------------
	# plot confidence or prediction bands
	if (ci_band) {
		# plot lme line and 95% confidence bands
		p = p + geom_ribbon(aes(x = get(x_var),y = get(y_var),ymin=get(y_var)-2*SE,ymax=get(y_var)+2*SE,
			colour = get(subgrp_var), fill = get(subgrp_var)),data = newdat, alpha = band_alpha) +
			geom_line(aes(x = get(x_var), y = get(y_var), colour = get(subgrp_var)),data = newdat)
	} else if (pi_band)	{
		# plot lme line and 95% prediction bands
		p = p + geom_ribbon(aes(x = get(x_var),y = get(y_var),ymin=get(y_var)-2*SE2,ymax=get(y_var)+2*SE2,
			colour = get(subgrp_var), fill = get(subgrp_var)),data = newdat, alpha = band_alpha) +
			geom_line(aes(x = get(x_var), y = get(y_var), colour = get(subgrp_var)),data = newdat)
	}# if (ci_band)

	# change y limits if necessary
	if (!is.null(yLimits)) {
		p = p + scale_y_continuous(limits = yLimits)
	}

	p = p + theme(legend.title=element_blank()) + xlab(xLabel) + ylab(yLabel) +
			scale_x_continuous(limits = xLimits)
	#------------------------------------------------------------------------------
	# save plot
	if (!is.null(fname2save)) {
		ggsave(filename = fname2save)
	}# if

	#------------------------------------------------------------------------------
	# standardize data
	if (standardize){
	  df2use[,y_var] = (df2use[,y_var] - mean(df2use[,y_var], na.rm=TRUE))/sd(df2use[,y_var], na.rm = TRUE)
	  df2use[,x_var] = (df2use[,x_var] - mean(df2use[,x_var], na.rm = TRUE))/sd(df2use[,x_var], na.rm = TRUE)

	  # re-run model on mean centered data
	  m_allgrps <- eval(substitute(lme(fixed = fx_form, random = rx_form, data = df2use,
	                                   na.action = na.omit, control=ctrl, method = ml_method),
	                               list(fx_form = fx_form, rx_form = rx_form)))
	}
	#------------------------------------------------------------------------------

	#------------------------------------------------------------------------------
	# output information
	p
	result = list(p = p, lme_model = m_allgrps)
} # spaghettiPlot
#==============================================================================



#==============================================================================
# function identify ET subgroup
getETsubgrp <- function(tidy_df, full_df, Dx) {
	# DESCRIPTION
	# Will find the ET subgroup for each subject and append it to tidy_df
	#
	# INPUT
	# tidy_df = the reduced df for longitudinal modeling
	# full_df = the full LW df in its raw form
	#

	tmp_df = full_df[,c("subjectId","p2f2","pct_Fixation_Geometric")]
	if (Dx=="ASD") {
		uniqueSubs = unique(tidy_df$subjectId)
		for (i in 1:length(uniqueSubs)) {
			tidy_sub_mask = is.element(tidy_df$subjectId,uniqueSubs[i])
			full_sub_mask = is.element(tmp_df$subjectId,uniqueSubs[i])
			tmp_pct_geo = tmp_df[full_sub_mask,"pct_Fixation_Geometric"]
			if (tmp_pct_geo>=69) {
				# print(sprintf("Geo ASD %f",tmp_pct_geo))
				tidy_df$ETsubgrpDx[tidy_sub_mask] = "Geo ASD"
			} else if (tmp_pct_geo<69) {
				# print(sprintf("nonGeo ASD %f",tmp_pct_geo))
				tidy_df$ETsubgrpDx[tidy_sub_mask] = "nonGeo ASD"
			} # if (tmp_pct_geo>=0.69)
		}# for (i in 1:length(uniqueSubs))
	} else {
		tidy_df$ETsubgrpDx = Dx
	}# if (Dx=="ASD")
	# tidy_df$ETsubgrpDx = as.factor(tidy_df$ETsubgrpDx)
	return(tidy_df)
}# function getETsubgrp
#==============================================================================




#==============================================================================
# function identify ET subgroup
getETsubgrp2 <- function(tidy_df, full_df, Dx) {
	# DESCRIPTION
	# Will find the ET subgroup for each subject and append it to tidy_df
	#
	# INPUT
	# tidy_df = the reduced df for longitudinal modeling
	# full_df = the full LW df in its raw form
	#

	tmp_df = full_df[,c("subjectId","p2f2","pct_Fixation_Geometric")]
	if (Dx=="ASD") {
		uniqueSubs = unique(tidy_df$subjectId)
		for (i in 1:length(uniqueSubs)) {
			#print(uniqueSubs[i])
			#print(i)
			tidy_sub_mask = is.element(tidy_df$subjectId,uniqueSubs[i])
			full_sub_mask = is.element(tmp_df$subjectId,uniqueSubs[i])
			tmp_pct_geo = tmp_df[full_sub_mask,"pct_Fixation_Geometric"]
			if (sum(full_sub_mask)==0) {
				# print(sprintf("nonGeo ASD %f",tmp_pct_geo))
				tidy_df$ETsubgrpDx[tidy_sub_mask] = "ASD"
			} else if (tmp_pct_geo>=69) {
				#print(sprintf("Geo ASD %f",tmp_pct_geo))
				tidy_df$ETsubgrpDx[tidy_sub_mask] = "Geo ASD"
			} else if (tmp_pct_geo<69) {
				# print(sprintf("nonGeo ASD %f",tmp_pct_geo))
				tidy_df$ETsubgrpDx[tidy_sub_mask] = "nonGeo ASD"
			} # if (tmp_pct_geo>=0.69)
		}# for (i in 1:length(uniqueSubs))
	} else {
		tidy_df$ETsubgrpDx = Dx
	}# if (Dx=="ASD")
	# tidy_df$ETsubgrpDx = as.factor(tidy_df$ETsubgrpDx)
	return(tidy_df)
}# function getETsubgrp2
