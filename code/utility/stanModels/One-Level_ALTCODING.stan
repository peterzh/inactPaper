// NOT COMPLETE, GETTING WEIRD BUG ON LINE 43

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
	real BL[numSessions];
	real BR[numSessions];
	real<lower=0> SL[numSessions];
	real<lower=0> SR[numSessions];
	real<lower=0> N_EXP[numSessions];

	cov_matrix[5] Sigma_sess; //covariance matrix for session variation
	cov_matrix[5] Sigma_subj; //covariance matrix for session variation
	
	vector[5] b_global;
	vector[5] b_sess[numSessions];
	vector[5] b_subj[numSubjects];
	

	}
transformed parameters {
	vector[3] logOdds[numTrials]; // ln pL/pNG, ln pR/pNG, ln pNG/pNG
	//compute (non)linear model

	for (n in 1:numTrials)
	{
		logOdds[n][1] = BL[sessionID[n]] + SL[sessionID[n]]*contrastLeft[n]^N_EXP[sessionID[n]];
		logOdds[n][2] = BR[sessionID[n]] + SR[sessionID[n]]*contrastRight[n]^N_EXP[sessionID[n]];
		logOdds[n][3] = 0;
	}
}
model {
	
	//hyperpriors on global parameters
	b_global[1:2] ~ normal(0, 2);
	b_global[3:4] ~ normal(5, 2);
	b_global[5] ~ normal(0.5,0.25);
	
	//prior on the parameters for each subject
	for (subj in 1:numSubjects) {
		b_subj[subj] ~ multi_normal(b_global, Sigma_subj);
	}
	

	//prior on the parameters for each session
	for (sess in 1:numSessions) {
		b_sess[sess][1] = BL[sess];
		b_sess[sess][2] = BR[sess];
		b_sess[sess][3] = SL[sess];
		b_sess[sess][4] = SR[sess];
		b_sess[sess][5] = N_EXP[sess];
		b_sess[sess] ~ multi_normal( b_subj[subjID_session[sess]], Sigma_sess);
	}
	
	
	for (n in 1:numTrials) {
		choice[n] ~ categorical_logit( logOdds[n] );
	}
}
generated quantities {
	vector[numTrials] log_lik;
	
	//write loglik
	for (n in 1:numTrials){
		log_lik[n] = categorical_logit_lpmf(choice[n] | logOdds[n] );
	} 
}