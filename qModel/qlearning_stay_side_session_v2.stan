data {
        int NT; //  MAX length of trials per session
        int NS; // max number of sessions per animal
        int NA; // number of animals
        int NS_all[NA]; // number of sessions per subject
        int NT_all[NA, NS]; // trial lengths per session per animal
        int NS_tot; // total number of sessions
        int r[NA,NS,NT]; // reward
        int c[NA,NS,NT]; // choice
        int subj_idx[NS_tot]; // NS_tot x 1 vector of subject IDs
        int sess_idx[NS_tot]; // NS_tot x 1 vector of session IDs (for each subject) 
        int start_idx[NA]; // NA x 1 vector of session indices for first session of the animal 
        }


parameters {

// Population level 
		// inverse temperature
        real betam;
        real<lower=0> betasd;
        real<lower=0> beta_sd_m; // mean and sd of the subject level sd
        real<lower=0> beta_sd_sd;
        
		// learning rate 
        real<lower=0> alpham;
        real<lower=0> alphasd;
        real<lower=0> alpha_sd_m;
        real<lower=0> alpha_sd_sd;
		
		// perseveration 
        real staym;
        real<lower=0> staysd;
		real<lower=0> stay_sd_m;
		real<lower=0> stay_sd_sd;
        
        // side bias
        real sidem;
        real<lower=0> sidesd;
		real<lower=0> side_sd_m;
		real<lower=0> side_sd_sd;
		
/// for each mouse
		real beta_mice_m[NA];
		real<lower=0>beta_mice_sd[NA];
		
		real alpha_mice_m[NA];
		real<lower=0>alpha_mice_sd[NA];
		
		real stay_mice_m[NA];
		real<lower=0>stay_mice_sd[NA];
		
		real side_mice_m[NA];
		real<lower=0>side_mice_sd[NA];
		
		
// for each mouse's session
		real betas[NS_tot];
		real<lower=0> alphas[NS_tot];
		real stays[NS_tot];
		real sides[NS_tot]; 
		}


model {

        betam ~ normal(0,2);
        alpham ~ normal(.5,1);
        staym ~ normal(0,1);
        sidem ~ normal(0,2);

        betasd ~ normal(0,2);
        alphasd ~ normal(0,2);
        staysd ~ normal(0,2);
        sidesd ~ normal(0,2);
        
        
        beta_sd_m ~ normal(0,2);
        alpha_sd_m ~ normal(0,2);
        stay_sd_m ~ normal(0,2);
        side_sd_m ~ normal(0,2);


        beta_sd_sd ~ normal(0,2);
        alpha_sd_sd ~ normal(0,2);
        stay_sd_sd ~ normal(0,2);
        side_sd_sd ~ normal(0,2);

			
        for (s in 1:NS_tot){
            real alpha;
            real q[2];
            real pc;
            	
            if (s == start_idx[subj_idx[s]]){
            	beta_mice_m[subj_idx[s]]  ~ normal(betam,betasd);
            	beta_mice_sd[subj_idx[s]] ~ normal(beta_sd_m, beta_sd_sd);
            	
            	alpha_mice_m[subj_idx[s]] ~ normal(alpham,alphasd);
				alpha_mice_sd[subj_idx[s]] ~ normal(alpha_sd_m, alpha_sd_sd); 
			
				stay_mice_m[subj_idx[s]] ~ normal(staym,staysd);
				stay_mice_sd[subj_idx[s]] ~ normal(stay_sd_m, stay_sd_sd); 
			
				side_mice_m[subj_idx[s]] ~ normal(sidem,sidesd);
				side_mice_sd[subj_idx[s]] ~ normal(side_sd_m, side_sd_sd);
            	
            	}
            	
            
        
        	
                
                
            betas[s] ~ normal(beta_mice_m[subj_idx[s]], beta_mice_sd[subj_idx[s]]);
            alphas[s] ~ normal(alpha_mice_m[subj_idx[s]], alpha_mice_sd[subj_idx[s]]);
            stays[s] ~ normal(stay_mice_m[subj_idx[s]], stay_mice_sd[subj_idx[s]]);
            sides[s] ~ normal(side_mice_m[subj_idx[s]], side_mice_sd[subj_idx[s]]);

            alpha <- Phi_approx(alphas[s]);


            for (i in 1:2) {
                    q[i] = 0;
                }


            for (t in 1:NT_all[subj_idx[s],sess_idx[s]]) {
                    if (t > 1) {
                            pc = c[subj_idx[s], sess_idx[s],(t - 1)] * 2 - 1;
                    } else {
                            pc = 0;
                    }

                    c[subj_idx[s],sess_idx[s],t] ~ bernoulli_logit(sides[s] + betas[s]  * (q[2] - q[1]) +  stays[s] * pc);


                    q[c[subj_idx[s],sess_idx[s],t]+1] = (1-alpha) * q[c[subj_idx[s],sess_idx[s],t]+1] + alpha * r[subj_idx[s],sess_idx[s],t]; // alpha taken out for rescale (added back to be \\  											//consistent with intercept_last

            }
	  }
} 
