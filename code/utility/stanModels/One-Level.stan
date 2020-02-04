data {
	int<lower=1> numTrials;
	int<lower=1> numSessions;
	int<lower=1,upper=numSessions> sessionID[numTrials];
	vector<lower=0,upper=1>[numTrials] contrastLeft;
	vector<lower=0,upper=1>[numTrials] contrastRight;
	int<lower=1,upper=3> choice[numTrials]; // 1=Left, 2=Right, 3=NoGo
}
parameters {
	//global parameters
	vector[2] bias; 
	vector<lower=0>[2] sens;
	real<lower=0> n_exp;
	
	//per session deviations 
	vector<lower=0>[2] sd_sess;
	matrix[2,numSessions] z_sess; //standard normal variable used to draw from the covariance matrix
	cholesky_factor_corr[2] rho_sess; //correlations of deviations
}
transformed parameters {
	vector[3] logOdds[numTrials]; // ln pL/pNG, ln pR/pNG, ln pNG/pNG
	matrix[2,numSessions] b_sess;
	
	//draw samples of sess and subj deviations, according to the covariance structure in rho_ & sd_
	b_sess = diag_pre_multiply(sd_sess, rho_sess) * z_sess;

	{
		//temp variables
		real BL;
		real BR;
		real SL;
		real SR;
		
		//compute (non)linear model
		for (n in 1:numTrials)
		{
			BL = bias[1] + b_sess[1,sessionID[n]];
			BR = bias[2] + b_sess[2,sessionID[n]];
			SL = sens[1];
			SR = sens[2];

			logOdds[n][1] = BL + SL*contrastLeft[n]^n_exp;		
			logOdds[n][2] = BR + SR*contrastRight[n]^n_exp;
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
	
	//prior on the variation of the per-session deviations
	sd_sess ~ cauchy(0,1);
	
	//prior on the cholesky factor of the covariance matrix
	rho_sess ~ lkj_corr_cholesky(2.0); //penalises extreme correlations between the deviations
	
	for (n in 1:numTrials) {
		choice[n] ~ categorical_logit( logOdds[n] );
	}
}
generated quantities {
	corr_matrix[2] corr_sess;
	vector[numTrials] log_lik;

	//write correlation matrix
	corr_sess = rho_sess * rho_sess';

	//write loglik
	for (n in 1:numTrials){
		log_lik[n] = categorical_logit_lpmf(choice[n] | logOdds[n] );
	} 
}