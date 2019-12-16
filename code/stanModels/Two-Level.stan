data {
	int<lower=1> numTrials;
	int<lower=1> numSessions;
	int<lower=1> numSubjects;
	int<lower=1,upper=numSessions> sessionID[numTrials];
	int<lower=1,upper=numSubjects> subjectID[numTrials];
	vector<lower=0,upper=1>[numTrials] contrastLeft;
	vector<lower=0,upper=1>[numTrials] contrastRight;
	int<lower=1,upper=3> choice[numTrials]; // 1=Left, 2=Right, 3=NoGo
	int<lower=1,upper=numSubjects> subjID_session[numSessions]; //the subject ID for each session in the dataset
}
parameters {
	//global parameters
	vector[2] bias; 
	vector<lower=0>[2] sens;
	real<lower=0> n_exp;
	
	//per session deviations 
	vector<lower=0>[4] sd_sess;
	matrix[4,numSessions] z_sess; //standard normal variable used to draw from the covariance matrix
	cholesky_factor_corr[4] rho_sess; //correlations of deviations
	
	//per subject deviations
	vector<lower=0>[4] sd_subj;
	matrix[4,numSubjects] z_subj; 
	cholesky_factor_corr[4] rho_subj; 
	}
transformed parameters {
	vector[3] logOdds[numTrials]; // ln pL/pNG, ln pR/pNG, ln pNG/pNG
	matrix[4,numSessions] b_sess;
	matrix[4,numSubjects] b_subj;
	
	//draw samples of sess and subj deviations, according to the covariance structure in rho_ & sd_
	b_sess = diag_pre_multiply(sd_sess, rho_sess) * z_sess;
	b_subj = diag_pre_multiply(sd_subj, rho_subj) * z_subj;

	{
		//temp variables
		real BL;
		real BR;
		real SL;
		real SR;
		real N_EXP;
		
		//compute (non)linear model
		for (n in 1:numTrials)
		{
			BL = bias[1] + b_sess[1,sessionID[n]] + b_subj[1,subjectID[n]];
			BR = bias[2] + b_sess[2,sessionID[n]] + b_subj[2,subjectID[n]];
			SL = sens[1] + b_sess[3,sessionID[n]] + b_subj[3,subjectID[n]];
			SR = sens[2] + b_sess[4,sessionID[n]] + b_subj[4,subjectID[n]];
			N_EXP = n_exp;

			logOdds[n][1] = BL + SL*contrastLeft[n]^N_EXP;
			logOdds[n][2] = BR + SR*contrastRight[n]^N_EXP;
			logOdds[n][3] = 0;
		}
	}
}
model {
	//priors on global parameters
	bias ~ normal(0, 2);
	sens ~ normal(5, 2);
	n_exp ~ normal(0.5,0.25);
	
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
	corr_matrix[4] corr_sess;
	corr_matrix[4] corr_subj;
	vector[numTrials] log_lik;

	
	//write correlation matrix
	corr_sess = rho_sess * rho_sess';
	corr_subj = rho_subj * rho_subj';

	//write loglik
	for (n in 1:numTrials){
		log_lik[n] = categorical_logit_lpmf(choice[n] | logOdds[n] );
	} 
}