data {
	int<lower=1> numTrials;
	vector<lower=0,upper=1>[numTrials] contrastLeft;
	vector<lower=0,upper=1>[numTrials] contrastRight;
	int<lower=1,upper=3> choice[numTrials]; // 1=Left, 2=Right, 3=NoGo
}
parameters {
	//global parameters
	vector[2] bias; 
	vector<lower=0>[2] sens;
	real<lower=0, upper=1> n_exp;
}
transformed parameters {
	vector[3] logOdds[numTrials]; // ln pL/pNG, ln pR/pNG, ln pNG/pNG
	{
		//compute (non)linear model
		for (n in 1:numTrials)
		{
			logOdds[n][1] = bias[1] + sens[1]*contrastLeft[n]^n_exp;		
			logOdds[n][2] = bias[2] + sens[2]*contrastRight[n]^n_exp;
			logOdds[n][3] = 0;
		}
	}
}
model {
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