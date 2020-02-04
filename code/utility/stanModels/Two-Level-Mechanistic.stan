data {
	int<lower=1> numTrials;
	int<lower=1> numSessions;
	int<lower=1> numSubjects;
	int<lower=1,upper=numSessions> sessionID[numTrials];
	int<lower=1,upper=numSubjects> subjectID[numTrials];
	int<lower=1,upper=3> choice[numTrials]; // 1=Left, 2=Right, 3=NoGo
	int<lower=1,upper=numSubjects> subjID_session[numSessions]; //the subject ID for each session in the dataset
	vector[4] firing_rate[numTrials];
}
parameters {
	//global parameters
	row_vector[2 + 8] w;
	
	//per session deviations 
	vector<lower=0>[2 + 8] sd_sess;
	matrix[2 + 8,numSessions] z_sess; //standard normal variable used to draw from the covariance matrix
	cholesky_factor_corr[2 + 8] rho_sess; //correlations of deviations
	
	//per subject deviations
	vector<lower=0>[2 + 8] sd_subj;
	matrix[2 + 8,numSubjects] z_subj; 
	cholesky_factor_corr[2 + 8] rho_subj; 
	}
transformed parameters {
	vector[3] logOdds[numTrials]; // ln pL/pNG, ln pR/pNG, ln pNG/pNG
	matrix[2 + 8,numSessions] b_sess;
	matrix[2 + 8,numSubjects] b_subj;
	
	//draw samples of sess and subj deviations, according to the covariance structure in rho_ & sd_
	b_sess = diag_pre_multiply(sd_sess, rho_sess) * z_sess;
	b_subj = diag_pre_multiply(sd_subj, rho_subj) * z_subj;

	{
			//temp variables
		real BL;
		real BR;
		row_vector[4] WL;
		row_vector[4] WR;
		
		for (n in 1:numTrials)
		{
			BL = w[1] + b_sess[1,sessionID[n]] + b_subj[1,subjectID[n]];
			BR = w[2] + b_sess[2,sessionID[n]] + b_subj[2,subjectID[n]];
			
			for (p in 1:4)
			{
				WL[p] = w[ 2 + 2*(p-1) + 1 ] + b_sess[2 + 2*(p-1) + 1,sessionID[n]] + b_subj[2 + 2*(p-1) + 1,subjectID[n]];
				WR[p] = w[ 2 + 2*(p-1) + 2 ] + b_sess[2 + 2*(p-1) + 2,sessionID[n]] + b_subj[2 + 2*(p-1) + 2,subjectID[n]];
			}
			
			logOdds[n][1] = BL + WL*firing_rate[n];
			logOdds[n][2] = BR + WR*firing_rate[n];
			logOdds[n][3] = 0;
		}
	}
}
model {
	//priors on global parameters
	w ~ normal(0, 4);
	
	//make z_std_normal be standard normal (non centred parameterisation)
	to_vector(z_sess) ~ normal(0, 1);	
	to_vector(z_subj) ~ normal(0, 1);	
	
	//prior on the variation of the per-session deviations
	sd_sess ~ cauchy(0,1);
	sd_subj ~ cauchy(0,1);
	
	//prior on the cholesky factor of the covariance matrix
	rho_sess ~ lkj_corr_cholesky(2.0); //penalises extreme correlations between the deviations
	rho_subj ~ lkj_corr_cholesky(2.0); //penalises extreme correlations between the deviations
	
	for (n in 1:numTrials) {
		choice[n] ~ categorical_logit( logOdds[n] );
	}
}
generated quantities {
	corr_matrix[2 + 8] corr_sess;
	corr_matrix[2 + 8] corr_subj;
	vector[numTrials] log_lik;
	
	//write correlation matrix
	corr_sess = rho_sess * rho_sess';
	corr_subj = rho_subj * rho_subj';

	//write loglik
	for (n in 1:numTrials){
		log_lik[n] = categorical_logit_lpmf(choice[n] | logOdds[n] );
	} 
}