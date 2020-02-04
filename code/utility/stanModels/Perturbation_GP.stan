//expanded perturbation model to account for the effect of pulsed inactivation
//at different times of the trial. Uses gaussian process to permit nearby times to have similar effects
data {
	int<lower=1> numTrials;
	vector<lower=0,upper=1>[numTrials] contrastLeft;
	vector<lower=0,upper=1>[numTrials] contrastRight;
	int<lower=1,upper=3> choice[numTrials]; // 1=Left, 2=Right, 3=NoGo
	int<lower=0,upper=1> perturbation[numTrials]; //yes-no perturbation
	real perturbationTime[numTrials]; // time when the perturbation happened in the trial
}
transformed data {
	int numPerturbationTrials;
	numPerturbationTrials = sum(perturbation);
}
parameters {
	//global parameters
	vector[2] bias; 
	vector<lower=0>[2] sens;
	real<lower=0> n_exp;

	//parameter perturbations by the laser
	real delta;
	vector[numPerturbationTrials] deltaT;
	
	//gaussian process parameters
	real<lower=0> eta_sq;
	real<lower=0> inv_rho_sq;
	real<lower=0> sigma_sq;
}
transformed parameters {
	vector[3] logOdds[numTrials]; // ln pL/pNG, ln pR/pNG, ln pNG/pNG

	//gp parameters
	real<lower=0> rho_sq;
	matrix[numPerturbationTrials,numPerturbationTrials] L; //cholesky factor
	cov_matrix[numPerturbationTrials] Sigma;
	vector[numPerturbationTrials] mu;
	rho_sq = inv(inv_rho_sq);
	for (i in 1:numPerturbationTrials) mu[i] = delta;
	
	//gp covariance definition
	for (i in 1:(numPerturbationTrials-1)) {
		for (j in (i+1):numPerturbationTrials) {
			Sigma[i, j] = eta_sq * exp(-rho_sq * pow(perturbationTime[i] - perturbationTime[j],2));
			Sigma[j, i] = Sigma[i, j];
		}
	}
	for (k in 1:numPerturbationTrials) {
		Sigma[k, k] = eta_sq + sigma_sq; // + jitter
	}

	L = cholesky_decompose(Sigma);

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
			SL = sens[1] ;
			SR = sens[2] ;
			N = n_exp;
			
			//bias perturbation, differs according to the TIME of the perturbation
			//but nearby times affect behaviour similarly, so we want model that similarity by a Gaussian Process
			if (perturbation[n]>0) 
			{
				//add inactivation effect at each time
				BL += deltaT[n];
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
	delta ~ normal(0,2);
	
	//gaussian process on perturbationTime effect
	eta_sq ~ cauchy(0, 5);
	inv_rho_sq ~ cauchy(0, 5);
	sigma_sq ~ cauchy(0, 5);
	
	deltaT ~ multi_normal_cholesky(mu, L);
	
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