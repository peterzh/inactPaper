data {
	int<lower=1> numTrials;
	int<lower=1> numSessions;
	int<lower=1> numSubjects;
	int<lower=1,upper=numSessions> sessionID[numTrials];
	int<lower=1,upper=numSubjects> subjectID[numTrials];
	vector<lower=0,upper=1>[numTrials] contrastLeft;
	vector<lower=0,upper=1>[numTrials] contrastRight;
	int<lower=1,upper=3> choice[numTrials]; // 1=Left, 2=Right, 3=NoGo
}
parameters {
	//global parameters
	real bias_GNG;
	real bias_LR;
	real<lower=0> sens_L;
	real<lower=0> sens_R;
	
}
transformed parameters {
	vector[3] pHat[numTrials];

	{
		//temp variables
		real C[2];
		
		real Z_GNG;
		real Z_LR;
		
		real P_GO;
		real P_LGO;
		real P_RGO;
		
		//compute model
		for (n in 1:numTrials)
		{
			C[1] = sens_L*contrastLeft[n];
			C[2] = sens_R*contrastRight[n];
			
			Z_GNG = bias_GNG + max(C);
			Z_LR = bias_LR + sens_L*contrastLeft[n] - sens_R*contrastRight[n];
			
			P_GO = exp(Z_GNG)/( 1 + exp(Z_GNG));
			P_LGO = exp(Z_LR)/( 1 + exp(Z_LR));
			P_RGO = 1 - P_LGO;
			
			pHat[n][1] = P_LGO * P_GO;
			pHat[n][2] = P_RGO * P_GO;
			pHat[n][3] = 1 - P_GO;
		}
	}
}
model {
	//priors on global parameters
	bias_GNG ~ normal(0, 3);
	bias_LR ~ normal(0, 3);
	sens_L ~ normal(0, 10);
	sens_R ~ normal(0, 10);
	
	for (n in 1:numTrials) {
		choice[n] ~ categorical( pHat[n] );
	}
}
generated quantities {
	vector[numTrials] log_lik;

	//write loglik
	for (n in 1:numTrials){
		log_lik[n] = categorical_lpmf(choice[n] | pHat[n] );
	} 

}