data {
	int<lower=1> numTrials;
	int<lower=1> numSessions;
	int<lower=1> numSubjects;
	int<lower=1,upper=numSessions> sessionID[numTrials];
	int<lower=1,upper=numSubjects> subjectID[numTrials];
	vector<lower=0,upper=1>[numTrials] contrastLeft;
	vector<lower=0,upper=1>[numTrials] contrastRight;
	int<lower=1,upper=3> choice[numTrials]; // 1=Left, 2=Right, 3=NoGo
	int<lower=0> numPerturbations; //number of perturbation conditions
	int<lower=0,upper=numPerturbations> perturbation[numTrials];
	int<lower=1,upper=numSubjects> subjID_session[numSessions]; //the subject ID for each session in the dataset
}
parameters {
	//global parameters
	vector[2] bias; 
	vector<lower=0>[2] sens;
	real<lower=0> n_exp;
	
	//per session deviations 
	vector<lower=0>[5] sd_sess;
	matrix[5,numSessions] z_sess; //standard normal variable used to draw from the covariance matrix
	cholesky_factor_corr[5] rho_sess; //correlations of deviations
	
	//per subject deviations
	vector<lower=0>[5] sd_subj;
	matrix[5,numSubjects] z_subj; 
	cholesky_factor_corr[5] rho_subj; 
	
	//parameter perturbations by the laser
	real delta[4,numPerturbations];
}
transformed parameters {
	vector[3] pHat[numTrials];
	matrix[5,numSessions] b_sess;
	matrix[5,numSubjects] b_subj;
	
	//draw samples of sess and subj deviations, according to the covariance structure in rho_ & sd_
	b_sess = diag_pre_multiply(sd_sess, rho_sess) * z_sess;
	b_subj = diag_pre_multiply(sd_subj, rho_subj) * z_subj;

	{
		//temp variables
		real BGNG;
		real BLR;
		real SL;
		real SR;
		real N;
		real C[2];
		real Z_GNG;
		real Z_LR;
		real P_GO;
		real P_LGO;
		real P_RGO;
		
		//compute (non)linear model
		for (n in 1:numTrials)
		{
			BGNG = bias[1] + b_sess[1,sessionID[n]] + b_subj[1,subjectID[n]];
			BLR = bias[2] + b_sess[2,sessionID[n]] + b_subj[2,subjectID[n]];
			SL = sens[1] + b_sess[3,sessionID[n]] + b_subj[3,subjectID[n]];
			SR = sens[2] + b_sess[4,sessionID[n]] + b_subj[4,subjectID[n]];
			N = n_exp 	 + b_sess[5,sessionID[n]] + b_subj[5,subjectID[n]];
			
			if (perturbation[n]>0) 
			{
				BGNG += delta[1,perturbation[n]];
				BLR += delta[2,perturbation[n]];
				SL += delta[3,perturbation[n]];
				SR += delta[4,perturbation[n]];
			}

			C[1] = SL*contrastLeft[n]^N;
			C[2] = SR*contrastRight[n]^N;
			
			Z_GNG = BGNG + max(C);
			Z_LR = BLR + C[1] - C[2];
			//P_GO = exp(Z_GNG)/( 1 + exp(Z_GNG));
			//P_LGO = exp(Z_LR)/( 1 + exp(Z_LR));
			P_GO = inv_logit(Z_GNG);
			P_LGO = inv_logit(Z_LR);
			P_RGO = 1 - P_LGO;
			
			pHat[n][1] = P_LGO * P_GO;
			pHat[n][2] = P_RGO * P_GO;
			pHat[n][3] = 1 - P_GO;
		}
	}
}
model {
	//priors on global parameters
	bias ~ normal(0, 3);
	sens ~ normal(0, 5);
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
	
	//priors on perturbations
	to_array_1d(delta) ~ normal(0,2);
	
	for (n in 1:numTrials) {
		choice[n] ~ categorical( pHat[n] );
	}
}
generated quantities {
	corr_matrix[5] corr_sess;
	corr_matrix[5] corr_subj;
	vector[numTrials] log_lik;
	vector[3] pTest[numTestContrasts,numSessions];

	//write correlation matrix
	corr_sess = rho_sess * rho_sess';
	corr_subj = rho_subj * rho_subj';

	//write loglik
	for (n in 1:numTrials){
		log_lik[n] = categorical_lpmf(choice[n] | pHat[n] );
	} 
	
	{
		//temp variables
		real BGNG;
		real BLR;
		real SL;
		real SR;
		real N;
		real C[2];
		real Z_GNG;
		real Z_LR;
		real P_GO;
		real P_LGO;
		real P_RGO;
		
		//write posterior predictions for each session
		for (sess in 1:numSessions) {
			BGNG = bias[1] + b_sess[1,sess] + b_subj[1,subjID_session[sess]];
			BLR = bias[2] + b_sess[2,sess] + b_subj[2,subjID_session[sess]];
			SL = sens[1] + b_sess[3,sess] + b_subj[3,subjID_session[sess]];
			SR = sens[2] + b_sess[4,sess] + b_subj[4,subjID_session[sess]];
			N = n_exp   + b_sess[5,sess] + b_subj[5,subjID_session[sess]];
			
			for (c in 1:numTestContrasts) {
				
				
				C[1] = SL*testContrastLeft[c]^N;
				C[2] = SR*testContrastRight[c]^N;
				
				Z_GNG = BGNG + max(C);
				Z_LR = BLR + C[1] - C[2];
				P_GO = exp(Z_GNG)/( 1 + exp(Z_GNG));
				P_LGO = exp(Z_LR)/( 1 + exp(Z_LR));
				P_RGO = 1 - P_LGO;
				
				pTest[c,sess][1] = P_LGO * P_GO;
				pTest[c,sess][2] = P_RGO * P_GO;
				pTest[c,sess][3] = 1 - P_GO;
				
				
			}
		}
		
		
	}

}