/*
-Hierarchical behavioural model for Pip's 2AFC data
-Also models parameter perturbations induced by optogenetic inactivation
-Perturbations to the parameters are NOT included in the hierarchy, all sessions/subjects are assumed to have the same delta for each parameter.
*/
data {
	int<lower=1> numTrials;
	int<lower=1> numSessions;
	int<lower=1> numSubjects;
	int<lower=1,upper=numSessions> sessionID[numTrials];
	int<lower=1,upper=numSubjects> subjectID[numTrials];
	vector<lower=0,upper=1>[numTrials] contrastLeft;
	vector<lower=0,upper=1>[numTrials] contrastRight;
	int<lower=0,upper=1> audioLeft[numTrials];
	int<lower=0,upper=1> audioRight[numTrials];
	int<lower=0,upper=1> audioZero[numTrials];
	int<lower=0,upper=1> choice[numTrials]; // 0=Left, 1=Right
}
parameters {
	//global parameters
	real bias;
	real<lower=0> vis;
	real<lower=0> aud_left;
	real<lower=0> aud_right;
	real<lower=0,upper=1> vis_n_exp;
	
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
	vector[numTrials] log_pRpL;
	matrix[4,numSessions] b_sess;
	matrix[4,numSubjects] b_subj;
	
	//draw samples of sess and subj deviations, according to the covariance structure in rho_ & sd_
	b_sess = diag_pre_multiply(sd_sess, rho_sess) * z_sess;
	b_subj = diag_pre_multiply(sd_subj, rho_subj) * z_subj;

	{
		//temp variables
		real B;
		real VL;
		real VR;
		real AL;
		real AR;
		real A0;
		
		//compute (non)linear model
		for (n in 1:numTrials)
		{
			B  = bias 		+ b_sess[1,sessionID[n]] + b_subj[1,subjectID[n]];
			VL = vis 		+ b_sess[2,sessionID[n]] + b_subj[2,subjectID[n]];
			VR = VL;
			AL = aud_left 	+ b_sess[3,sessionID[n]] + b_subj[3,subjectID[n]];
			AR = aud_right 	+ b_sess[4,sessionID[n]] + b_subj[4,subjectID[n]];
			
			log_pRpL[n] = B + VR*(contrastRight[n]^vis_n_exp) + AR*audioRight[n] - VL*(contrastLeft[n]^vis_n_exp) - AL*audioLeft[n];
		}
	}
}
model {
	//priors on global parameters
	bias 		~ normal(0, 2);
	vis 		~ normal(4, 3);
	aud_left 	~ normal(2, 3);
	aud_right 	~ normal(2, 3);
	vis_n_exp	~ normal(0.6,0.3);
	
	//make z_std_normal be standard normal (non centred parameterisation)
	to_vector(z_sess) ~ normal(0, 1);	
	to_vector(z_subj) ~ normal(0, 1);	
	
	//prior on the variation of the per-session deviations
	sd_sess ~ cauchy(0,1);
	sd_subj ~ cauchy(0,1);
	
	//prior on the cholesky factor of the covariance matrix
	rho_sess ~ lkj_corr_cholesky(2.0); //penalises extreme correlations between the deviations
	rho_subj ~ lkj_corr_cholesky(2.0); //penalises extreme correlations between the deviations

	//likelihood function
	choice ~ bernoulli_logit( log_pRpL );
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
		log_lik[n] = bernoulli_logit_lpmf(choice[n] | log_pRpL[n] );
	} 
}