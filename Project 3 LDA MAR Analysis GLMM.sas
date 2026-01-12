data alzheimer25;
    set "/home/u64400254/sasuser.v94/alzheimer25.sas7bdat"; 
run;

proc means data=alzheimer25 noprint;
    var AGE BMI;
    output out=stats mean=mean_age mean_bmi std=sd_age sd_bmi;
run;

data alzheimer_long;
    if _n_ = 1 then set stats;   
    set alzheimer25;
    
    array cdrsb_arr[0:6] cdrsb0-cdrsb6;
    array bprs_arr[0:6] bprs0-bprs6;
    array abpet_arr[0:6] abpet0-abpet6;
    array taupet_arr[0:6] taupet0-taupet6;
    
    do TIME = 0 to 6;
        CDRSB  = cdrsb_arr[TIME];
        BPRS   = bprs_arr[TIME];
        ABPET  = abpet_arr[TIME];
        TAUPET = taupet_arr[TIME];

        if not missing(CDRSB) 
         and not missing(BPRS)
         and not missing(ABPET)
         and not missing(TAUPET) then do;

            if CDRSB <= 10 then CDRSB_CAT = 0;
            else CDRSB_CAT = 1;

            AGE_STD = (AGE - mean_age) / sd_age;
            BMI_STD = (BMI - mean_bmi) / sd_bmi;

        	TIMECLSS = put(TIME, 1.);

            output;
        end;
    end;

	keep PATID SEX AGE AGE_STD BMI BMI_STD ADL 
     	TIME TIMECLSS CDRSB CDRSB_CAT BPRS ABPET CDRSB0 TAUPET ABPET0 TAUPET0
     	TRIAL JOB WZC EDU INKOMEN;
run;

proc means data=ALZHEIMER_LONG noprint;
    var ABPET TAUPET; 
    output out=bio_stats mean=m_ab m_tau std=s_ab s_tau;
run;

data alzheimer_long_centered;
    if _n_=1 then set bio_stats;
    set alzheimer_long;
    
    ABPET_STD = (ABPET - m_ab) / s_ab;
    TAUPET_STD = (TAUPET - m_tau) / s_tau;
run;

/* ============================================================== */

/* GLMM - Generalized Linear Mixed-Effects Model */

proc glimmix data=alzheimer_long_centered method=QUAD(QPOINTS=10);
    class PATID SEX;
    
    model CDRSB_CAT(descending) = TIME AGE_STD BMI_STD SEX ADL ABPET_STD TAUPET_STD
	TIME*BMI_STD TIME*TAUPET_STD TIME*SEX TIME*AGE_STD TIME*ABPET_STD TIME*ADL
          / dist=binary link=logit solution;
          
    random intercept TIME/ subject=PATID TYPE=UN solution;
    
	ods output SolutionR=eb_predictions;
    output out=glmm_preds pred(ilink)=Predicted_Prob;
    
    title "Model 3: Generalized Linear Mixed-Effects Model";
run;

/* ============================================================== */

/* Multiple Imputation */

/* The missing data at each time point is predicted using 
regression models conditioned on previous observed values. */
proc mi data=alzheimer25 out=mi_wide nimpute=10 seed=3 noprint;
    class SEX; 
   	/* These gives warnings that indicate high multicollinearity or data sparsity (few patients) at later visits. 
   	SAS automatically drops redundant predictors to allow the imputation model to converge. */
    var SEX AGE BMI ADL 
        cdrsb0 taupet0 abpet0
        cdrsb1 taupet1 abpet1
        cdrsb2 taupet2 abpet2
        cdrsb3 taupet3 abpet3
        cdrsb4 taupet4 abpet4
        cdrsb5 taupet5 abpet5
        cdrsb6 taupet6 abpet6; 
	monotone reg;
run;

proc sort data=mi_wide; by PATID _Imputation_; run;

data mi_long_glmm;
    set mi_wide;
    
    if _n_=1 then set stats;      
    if _n_=1 then set bio_stats;  
    
    array cdrsb_arr[0:6] cdrsb0-cdrsb6;
    array taupet_arr[0:6] taupet0-taupet6;
    array abpet_arr[0:6] abpet0-abpet6;

    do TIME = 0 to 6;
        CDRSB = cdrsb_arr[TIME];
        TAUPET = taupet_arr[TIME];
        ABPET = abpet_arr[TIME];

        if CDRSB < 10 then CDRSB_CAT = 0; 
        else CDRSB_CAT = 1;

        BMI_STD    = (BMI - mean_bmi) / sd_bmi;
        TAUPET_STD = (TAUPET - m_tau) / s_tau;
        ABPET_STD  = (ABPET - m_ab) / s_ab;
        AGE_STD    = (AGE - mean_age) / sd_age;
        
        output;
    end;
    
    keep _Imputation_ PATID TIME CDRSB_CAT BMI_STD TAUPET_STD 
         AGE_STD SEX ADL ABPET_STD;
run;

proc sort data=mi_long_glmm; by _Imputation_ PATID; run;

proc glimmix data=mi_long_glmm method=RSPL; /*QUAD*/
    by _Imputation_; 
    
    class PATID; 
    
    model CDRSB_CAT(descending) = TIME AGE_STD BMI_STD SEX ADL ABPET_STD TAUPET_STD
	TIME*BMI_STD TIME*TAUPET_STD TIME*SEX TIME*AGE_STD TIME*ABPET_STD TIME*ADL
          / dist=binary link=logit solution;
          
    random intercept TIME / subject=PATID TYPE=UN;
    
    ods output ParameterEstimates=glmm_parms;
    
    output out=glmm_predictions 
           pred(ilink)=PredProb 
           lcl(ilink)=LowerCI 
           ucl(ilink)=UpperCI;
run;

data glmm_ready;
    set glmm_parms;
    length VarName $50;
    
    VarName = compress(translate(Effect, "_", "*"), " ");
    
    keep _Imputation_ VarName Estimate StdErr;
run;

proc sort data=glmm_ready; by _Imputation_; run;

proc transpose data=glmm_ready out=betas(drop=_Name_);
    by _Imputation_;
    id VarName;
    var Estimate;
run;

/* Standard Errors (SE) */
proc transpose data=glmm_ready out=errors_se(drop=_Name_) prefix=SE_;
    by _Imputation_;
    id VarName;
    var StdErr;
run;

data combined_final;
    merge betas errors_se;
    by _Imputation_;
run;

/* Variance Between: variance of the estimate - Variance Within: mean(SE_^2) */
proc mianalyze data=combined_final;

    modeleffects Intercept 
                 TIME AGE_STD BMI_STD SEX ADL ABPET_STD TAUPET_STD
                 TIME_BMI_STD TIME_TAUPET_STD TIME_SEX 
                 TIME_AGE_STD TIME_ABPET_STD TIME_ADL;

    stderr       SE_Intercept 
                 SE_TIME SE_AGE_STD SE_BMI_STD SE_SEX SE_ADL SE_ABPET_STD SE_TAUPET_STD
                 SE_TIME_BMI_STD SE_TIME_TAUPET_STD SE_TIME_SEX 
                 SE_TIME_AGE_STD SE_TIME_ABPET_STD SE_TIME_ADL;
                 
    title "GLMM Results after Multiple Imputation";
run;

/* ============================================================== */

/* GLMM Pooled Prediction */

proc sql;
    create table plot_data as
    select TIME, 
           mean(PredProb) as Mean_Prob,   
           mean(LowerCI) as Mean_Lower,   
           mean(UpperCI) as Mean_Upper    
    from glmm_predictions
    group by TIME;
quit;

title "GLMM Pooled Prediction";
proc sgplot data=plot_data;

    band x=TIME upper=Mean_Upper lower=Mean_Lower / 
         transparency=0.5 legendlabel="95% Confidence Interval";
    

    series x=TIME y=Mean_Prob / 
           lineattrs=(thickness=2 color=blue) legendlabel="Predicted Probability";
    
    yaxis label="Prob(CDRSB_CAT = 1)" min=0 max=1 grid;
    xaxis label="Year" grid;
run;

/* ============================================================== */

/* Fitted Variance */

/* Var = (p * (1-p)) */
data plot_binomial_fitted;
    set plot_data; 
    
    Fitted_Var_Binomial = Mean_Prob * (1 - Mean_Prob);
run;

title "GLMM with MI Fitted Variance Function";
title2 "Variance = p * (1-p)";

proc sgplot data=plot_binomial_fitted;

    series x=TIME y=Fitted_Var_Binomial / 
           lineattrs=(color=blue thickness=2) 
           legendlabel="Fitted Variance";

    yaxis label="Variance" grid min=0.2 max=0.25; 
    xaxis label="Years" grid values=(0 to 6 by 1);
run;

/* ============================================================== */

proc means data=alzheimer_long noprint nway;
    class TIME;
    var CDRSB_CAT;
    output out=empirical_stats var=Empirical_Var;
run;

data plot_comparison;
    merge plot_binomial_fitted empirical_stats;
    by TIME;
run;

title "Empirical vs Fitted Variance GLMM-MI";

proc sgplot data=plot_comparison;
    /* GLMM-MI Variance */
    series x=TIME y=Fitted_Var_Binomial / 
           lineattrs=(color=blue thickness=2) 
           legendlabel="Fitted Variance (Theoretical)";

    /* Empirical Variance */
    series x=TIME y=Empirical_Var / 
            lineattrs=(color=red thickness=2) 
            legendlabel="Empirical Variance (Observed)";

    yaxis label="Variance" grid; 
    xaxis label="Years" grid values=(0 to 6 by 1);
run;
