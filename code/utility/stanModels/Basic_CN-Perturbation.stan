data {
	int<lower=1> numTrials;
	vector<lower=0,upper=1>[numTrials] contrastLeft;
	vector<lower=0,upper=1>[numTrials] contrastRight;
	int<lower=1,upper=3> choice[numTrials]; // 1=Left, 2=Right, 3=NoGo
	int<lower=0> numPerturbations; //number of perturbation conditions
	int<lower=0,upper=numPerturbations> perturbation[numTrials];
}
parameters {
	//global parameters
	vector[2] bias; 
	vector<lower=0>[2] sens;
	real<lower=0> n_exp;

	//parameter perturbations by the laser
	real delta[4,numPerturbations];
}
transformed parameters {
	vector[3] logOdds[numTrials]; // ln pL/pNG, ln pR/pNG, ln pNG/pNG

	{
		//temp variables
		real BL;
		real BR;
		real SL;
		real SR;
		real N;
		
		//compute (non)linear model
		for (n in 1:numTrials)
		{
			BL = bias[1];
			BR = bias[2];
			SL = sens[1];
			SR = sens[2];
			N = n_exp;
			
			if (perturbation[n]>0) 
			{
				BL += delta[1,perturbation[n]];
				BR += delta[2,perturbation[n]];
				SL += delta[3,perturbation[n]];
				SR += delta[4,perturbation[n]];
			}

			logOdds[n][1] = BL + SL*contrastLeft[n]^N;		
			logOdds[n][2] = BR + SR*contrastRight[n]^N;
			logOdds[n][3] = 0;
		}
	}
}
model {
	//priors on global parameters
	bias ~ normal(0, 2);
	sens ~ normal(5, 2);
	n_exp ~ normal(0.5,0.25);
	
	//priors on perturbations
	to_array_1d(delta) ~ normal(0,2);
	
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