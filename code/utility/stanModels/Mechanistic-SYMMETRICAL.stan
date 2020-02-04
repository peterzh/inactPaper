data {
	int<lower=1> numTrials;
	int<lower=1> numSessions;
	int<lower=1> numSubjects;
	int<lower=1,upper=3> choice[numTrials]; // 1=Left, 2=Right, 3=NoGo
	vector[4] firing_rate[numTrials];
}
parameters {
	//global parameters
	row_vector[2 + 4] w;
	}
transformed parameters {
	vector[3] logOdds[numTrials]; // ln pL/pNG, ln pR/pNG, ln pNG/pNG
	
	{
			//temp variables
		real BL;
		real BR;
		row_vector[4] W;
		int WR_IDX[4];
		WR_IDX[1] = 2;
		WR_IDX[2] = 1;
		WR_IDX[3] = 4;
		WR_IDX[4] = 3;
				
		for (n in 1:numTrials)
		{
			BL = w[1];
			BR = w[2];
			
			for (p in 1:4)
			{
				W[p] = w[2+p];
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