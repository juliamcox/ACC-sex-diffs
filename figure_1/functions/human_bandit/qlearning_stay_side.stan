data {
        int NT; //  MAX length of trials
        int NS; // number of subjects
        int NT_all[NS]; // trial lengths per subject
        int r[NS,NT]; //reward
        int c[NS,NT]; //choice
}


parameters {

        real betam;
        real alpham;
        real staym;
        real sidem;


        real<lower=0> betasd;
        real<lower=0> alphasd;
        real<lower=0> staysd;
        real<lower=0> sidesd;



        real betas[NS];
        real alphas[NS];
        real stays[NS];
        real sides[NS];

}


model {

        betam ~ normal(0,2);
        alpham ~ normal(0,2);
        staym ~ normal(0,2);
        sidem ~ normal(0,2);

        betasd ~ normal(0,2);
        alphasd ~ normal(0,2);
        staysd ~ normal(0,2);
        sidesd ~ normal(0,2);



        for (s in 1:NS) {

        	real alpha;
                real pc;
                real q[2];

                betas[s] ~ normal(betam, betasd);
                alphas[s] ~ normal(alpham, alphasd);
                stays[s] ~ normal(staym, staysd);
                sides[s] ~ normal(sidem, sidesd);

                alpha <- Phi_approx(alphas[s]);


                for (i in 1:2) {
                        q[i] = .5;
                }


                for (t in 1:NT_all[s]) {
                        if (t > 1) {
                                pc = c[s, (t - 1)] * 2 - 1;
                        } else {
                                pc = 0;
                        }

                        c[s,t] ~ bernoulli_logit(sides[s] + betas[s]  * (q[2] - q[1]) +  stays[s] * pc);


                        q[c[s,t]+1] = (1-alpha) * q[c[s,t]+1] + alpha * r[s,t]; // alpha taken out for rescale (added back to be consistent with intercept_last

                }

        }

}

