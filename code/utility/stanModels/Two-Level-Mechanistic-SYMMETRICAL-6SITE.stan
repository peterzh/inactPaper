data {
	int<lower=1> numTrials;
	int<lower=1> numSessions;
	int<lower=1> numSubjects;
	int<lower=1,upper=numSessions> sessionID[numTrials];
	int<lower=1,upper=numSubjects> subjectID[numTrials];
	int<lower=1,upper=3> choice[numTrials]; // 1=Left, 2=Right, 3=NoGo
	int<lower=1,upper=numSubjects> subjID_session[numSessions]; //the subject ID for each session in the dataset
	vector[6] firing_rate[numTrials];
}
parameters {
	//global parameters
	row_vector[2 + 6] w;
	
	//per session deviations 
	vector<lower=0>[2 + 6] sd_sess;
	matrix[2 + 6,numSessions] z_sess; //standard normal variable used to draw from the covariance matrix
	cholesky_factor_corr[2 + 6] rho_sess; //correlations of deviations
	
	//per subject deviations
	vector<lower=0>[2 + 6] sd_subj;
	matrix[2 + 6,numSubjects] z_subj; 
	cholesky_factor_corr[2 + 6] rho_subj; 
	}
transformed parameters {
	vector[3] logOdds[numTrials]; // ln pL/pNG, ln pR/pNG, ln pNG/pNG
	matrix[2 + 6,numSessions] b_sess;
	matrix[2 + 6,numSubjects] b_subj;
	
	//draw samples of sess and subj deviations, according to the covariance structure in rho_ & sd_
	b_sess = diag_pre_multiply(sd_sess, rho_sess) * z_sess;
	b_subj = diag_pre_multiply(sd_subj, rho_subj) * z_subj;

	{
			//temp variables
		real BL;
		real BR;
		row_vector[6] W;
		int WR_IDX[6];
		WR_IDX[1] = 2;
		WR_IDX[2] = 1;
		WR_IDX[3] = 4;
		WR_IDX[4] = 3;
		WR_IDX[5] = 6;
		WR_IDX[6] = 5;
				
		for (n in 1:numTrials)
		{
			BL = w[1] + b_sess[1,sessionID[n]] + b_subj[1,subjectID[n]];
			BR = w[2] + b_sess[2,sessionID[n]] + b_subj[2,subjectID[n]];
			
			for (p in 1:6)
			{
				W[p] = w[2+p] + b_sess[2+p,sessionID[n]] + b_subj[2+p,subjectID[n]];
			}
			
			logOdds[n][1] = BL + W*firing_rate[n];
			logOdds[n][2] = BR + W[WR_IDX]*firing_rate[n];
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
	corr_matrix[2 + 6] corr_sess;
	corr_matrix[2 + 6] corr_subj;
	vector[numTrials] log_lik;
	
	//write correlation matrix
	corr_sess = rho_sess * rho_sess';
	corr_subj = rho_subj * rho_subj';

	//write loglik
	for (n in 1:numTrials){
		log_lik[n] = categorical_logit_lpmf(choice[n] | logOdds[n] );
	} 
}