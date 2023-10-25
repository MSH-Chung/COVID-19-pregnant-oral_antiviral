*** VARIABLES
* pseudo_key: pseudo id of patients
* group: 1 refers to nirmatrelvir/ritonavir users; and 0 refers to controls
* age: maternal age of patients (year)
* vaccine_full: COVID-19 vaccination status of patients (fully vaccinated vs not-fully vaccinated/ unvaccinated)
* vaccine_status_cate: COVID-19 vaccination status of patients (boosted vs non-boosted)

*** Study outcomes
* covid_hosp: COVID-19-related hospitalization outcome
* MMMI: MMMI outcome
* vag_bleed: vaginal bleeding outcome
* preg_ht: pregnancy-induced hypertension outcome
* preeclampsia: pre-eclampsia outcome
* eclampsia: eclampsia outcome
* HELLP: HELLP syndrome outcome
* preterm_birth: preterm_birth (28-day) outcome
* unspec_comp: unspecified complication of pregnancy outcome
* abx: infection requiring antibiotics outcome
* ICU: admission to ICU outcome
* HDU: referral to HDU outcome
* maternal_death: maternal death outcome
* c_section: caesarean section outcome
* preterm_birth_comb: preterm birth outcome
* stillbirth: stillbirth outcome
* nnd: neonatal death outcome

********************************************************************************
use "covid pregnant antiviral outpatient.dta", clear
keep pseudo_key group date_covid age age_gp gender vaccine_status charlson_index asthma cancer cvd lung obesity dm t1dm t2dm immuno_hist dm ht liver lung chd ckd cancer ym_covid_cate month_covid HELLP_bl GDM MMMI_bl preg_ht_preeclampsia_bl unspec_comp_bl covid_hist
compress
save "covid pregnant antiviral outpatient characteristics.dta", replace


**# ***** 1:10 PS matching
noi forvalues k = 1/7 {
	use "covid pregnant antiviral outpatient.dta", clear
	gen covid_infect = (date_covid >= mdy(9,1,2022))
	gen subgp_1 = 1
	gen subgp_2 = age <= 32
	gen subgp_3 = age > 32
	gen subgp_4 = vaccine_full == 1
	gen subgp_5 = vaccine_full == 0
	gen subgp_6 = vaccine_status_cate == 1
	gen subgp_7 = vaccine_status_cate == 2
	keep if subgp_`k' == 1
	local mi = "logit group age c.date_covid#i.covid_infect i.vaccine_status age##i.MMMI_bl i.immuno_hist i.HELLP_bl i.GDM i.MMMI_bl i.preg_ht_preeclampsia_bl"
	
	`mi'
	predict prob_group
	gen xb_group = ln(prob_group / (1-prob_group))

	rename group group_old
	noi di "matching 1:10"
	set seed 123456
	calipmatch, gen(group_10_match) case(group_old) max(10) calipermatch(xb_group) caliperwidth(.05)
	bysort pseudo_key (group_10_match): replace group_10_match = group_10_match[1]
	gen group_10 = group_old if group_10_match < .

	compress
	save "covid pregnant antiviral outpatient matched subgp_`k'.dta", replace
}
*
********************************************************************************
**# Table 1 - Baseline characteristics after 1:10 ps matching
* N/mean, %/SD
qui forvalues j = 1/1 {
	forvalues k = 1(-1)0 {
		use "covid pregnant antiviral outpatient matched subgp_`j'.dta", clear
		replace group_10 = . if date_death <= date_baseline
		replace name_dose_1_num = . if name_dose_1_num == 4
		gen covid_infect0 = 1 - covid_infect
		keep if group_10 == `k'
		noi di _newline "after"
		noisily di "group_10=`k'"

		noisily di "var" _col(25) "N" _col(35) "mean" _col(50) "sd"
		foreach var in age {
			sum `var', d
			scalar m1 = r(mean)
			scalar n1 = r(N)
			scalar sd1 = r(sd)
			noi di "`var'" _col(25) n1 _col(35) m1 _col(50) sd1
		}

		noisily di "var" _col(25) "N" _col(35) "%"
		foreach var in covid_infect0 covid_infect vaccine_full ///
			chd ckd liver cancer GDM cholestasis MMMI_bl vag_bleed_bl preg_ht_bl preeclampsia_bl eclampsia_bl HELLP_bl unspec_comp_bl {
			sum `var', d
			scalar m1 = r(mean)
			scalar n1 = r(N)
			scalar b = round(n1*m1)
			noisily di "`var'" _col(25) b _col(35) m1*100
		}
		foreach var in age_gp vaccine_status_cate {
	    	noi tab `var'
		}
	}
}
*	
* SMD
qui	forvalues j = 1/1 {
	use "covid pregnant antiviral outpatient matched subgp_`j'.dta", clear
	replace group_10 = . if date_death <= date_baseline
	replace name_dose_1_num = . if name_dose_1_num == 4
	replace age_gp = age > 32
	noi di _newline "var" _col(25) "SMD"
***continuous
	foreach var in age {
		stddiff `var', by(group_10) abs cohensd
		noisily di "`var'" _col(25) abs(r(stddiff)[1,1])
	}
	noi di
***categorical
* binary
	foreach var in covid_infect vaccine_full ///
	chd ckd liver cancer GDM cholestasis MMMI_bl vag_bleed_bl preg_ht_bl preeclampsia_bl eclampsia_bl HELLP_bl unspec_comp_bl ///
	age_gp vaccine_status_cate {
		tab `var' group_10
		if r(r) > 1 {
			stddiff i.`var', by(group_10) abs cohensd
			noisily di "`var'" _col(25) abs(r(stddiff)[1,1])
		}
		else {
			noisily di "`var'" _col(25) .
		}
	}
	noisily di
}
*

********************************************************************************
***Creates directory
capture mkdir "matched_10"
forvalues k = 1/7 {
	capture mkdir "matched_10/subgp`k'"
	capture mkdir "matched_10/subgp`k'/prepare"
	capture mkdir "matched_10/subgp`k'/cloned"
	capture mkdir "matched_10/subgp`k'/split"
	capture mkdir "matched_10/subgp`k'/Treatment"
	capture mkdir "matched_10/subgp`k'/NoTreatment"
	capture mkdir "matched_10/subgp`k'/emulated"
	capture mkdir "matched_10/subgp`k'/weight"
	capture mkdir "matched_10/subgp`k'/KM ipcw"
}
*

**# ***** Target trial emulation
***** Prepare dataset for analysis
cls
* COVID-19-related hospitalization, MMMI, and individual MMMI components
* Follow-up of 28 days
qui foreach event in covid_hosp MMMI vag_bleed preg_ht preeclampsia eclampsia HELLP preterm_birth unspec_comp abx ICU HDU maternal_death {
forvalues k = 1/7 {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_`k'" _col(40) "bs`j'" _col(55) "setup" _col(70) "matched_10"
	}
	* Setup for trial emulation
	use "covid pregnant antiviral outpatient matched subgp_`k'.dta", clear
	rename group_10 group
	keep if group < .
	
	* events
	merge 1:1 pseudo_key using "date_covid_hosp", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_MMMI", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_ICU_HDU", keep(1 3) nogen
	
	if `j' > 0 {
		set seed `j'
		bsample
	}
	replace date_paxlovid = . if date_paxlovid >= date_admission
	* organize / rename / generate variables
	gen date_last_fu = min(date_death, mdy(02,12,2023), date_baseline + 28)
	gen date_event = date_`event'
	gen event = inrange(date_event, date_baseline, date_last_fu)
	gen fup_obs = min(date_event-date_baseline, 28) if event == 1
	replace fup_obs = min(date_last_fu-date_baseline, 28) if event == 0
	tab fup_obs
	gen time_to_treatment = date_paxlovid - date_baseline
	tab time_to_treatment
	gen treatment = inrange(date_paxlovid - date_baseline, 0, 5)
	* keep necessary variables
	keep pseudo_key fup_obs event time_to_treatment treatment
	gen bs = `j'
	compress
	save "matched_10/subgp`k'/prepare/prepare `event' subgp_`k' bs`j'.dta", replace
}
}
}
*
* Caesarean section, preterm birth, stillbirth, and neonatal death
* Follow-up not restricted to 28 days
qui foreach event in c_section preterm_birth_comb stillbirth nnd {
forvalues k = 1/7 {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_`k'" _col(40) "bs`j'" _col(55) "setup" _col(70) "matched_10"
	}
	* Setup for trial emulation
	use "covid pregnant antiviral outpatient matched subgp_`k'.dta", clear
	rename group_10 group
	keep if group < .
	keep if subgp_`k' == 1
	
	* events
	merge 1:1 pseudo_key using "date_c_section", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_preterm_birth", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_stillbirth", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_neonatal_death", keep(1 3) nogen
	
	if `j' > 0 {
		set seed `j'
		bsample
	}
	replace date_paxlovid = . if date_paxlovid >= date_admission
	* organize / rename / generate variables
	gen date_last_fu = min(date_death, mdy(02,12,2023))
	gen date_event = date_`event'
	gen event = inrange(date_event, date_baseline, date_last_fu)
	gen fup_obs = date_event-date_baseline if event == 1
	replace fup_obs = date_last_fu-date_baseline if event == 0
	tab fup_obs
	gen time_to_treatment = date_paxlovid - date_baseline
	tab time_to_treatment
	gen treatment = inrange(date_paxlovid - date_baseline, 0, 5)
	* keep necessary variables
	keep pseudo_key fup_obs event time_to_treatment treatment
	gen bs = `j'
	compress
	save "matched_10/subgp`k'/prepare/prepare `event' subgp_`k' bs`j'.dta", replace
}
}
}
*
cls
***** Bootstrap
*** Cloning & censoring
qui foreach event in covid_hosp MMMI vag_bleed preg_ht preeclampsia eclampsia HELLP preterm_birth unspec_comp abx ICU HDU maternal_death c_section preterm_birth_comb stillbirth nnd {
forvalues k = 1/7 {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_`k'" _col(40) "bs`j'" _col(55) "cloning" _col(70) "matched_10" 
	}
	* Prepare dataset for analysis
	use "matched_10/subgp`k'/prepare/prepare `event' subgp_`k' bs`j'.dta", replace
	stset fup_obs, failure(event)

	* Arm A: no treatment within 5 days (control: non-exposed group)
	gen outcomeA = _d // _d = `event'
	gen fupA = _t // _t = follow up time

	/// if the patient received treatment within 5 days:
	/// 1. no event outcome, since the patient survived till censoring (treatment)
	replace outcomeA = 0 if treatment==1 & time_to_treatment <=5 
	/// 2. follow up is censored at treatment
	replace fupA = time_to_treatment if treatment==1 & time_to_treatment <=5

	* Arm B: treatment within 5 days (treated: exposed group)
	gen outcomeB = _d 
	gen fupB = _t 

	/// if the patient survived the first 5 days and did not receive treatment within 5 days:
	/// 1. no event outcome if the patient survived the first 5 days
	replace outcomeB = 0 if (treatment==0 & _t>5) | (treatment==1 & time_to_treatment >5 & time_to_treatment !=.)
	/// 2. follow up is censored at 5 days
	replace fupB = 5 if (treatment==0 & _t>5) | (treatment==1 & time_to_treatment >5 & time_to_treatment != .)

	** append clones 
	preserve
		drop outcomeB fupB
		rename outcomeA outcome
		rename fupA fup
		gen arm = "NoTreatment"
		tempfile a
		save "`a'", replace
	restore
		drop outcomeA fupA
		rename outcomeB outcome
		rename fupB fup
		gen arm = "Treatment"	
		cap append using "`a'"

	// Weight models

	sort _all
	gen NewID = _n

	** add 1 day to 0-survivors
	replace fup= 1 if fup==0

	** Weight model: define survival time and event indicator	
	* treatment Arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	gen wm_fup = time_to_treatment if arm == "Treatment" & time_to_treatment<=5 & time_to_treatment!=. & treatment == 1 
	gen wm_outcome = 0 if arm == "Treatment" & time_to_treatment<=5 & time_to_treatment!=. & treatment == 1 

	** Case 2: they deviate at 5 days
	replace wm_fup = 5 if arm == "Treatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1))
	replace wm_outcome = 1 if arm == "Treatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "Treatment" & treatment == 0 & fup < 5
	replace wm_outcome = 0 if arm == "Treatment" & treatment == 0 & fup < 5
	** add 1 days to 0-survivors
	replace wm_fup= 1 if arm == "Treatment" & wm_fup==0 

	* No treatment Arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	replace wm_fup = time_to_treatment if arm == "NoTreatment" & time_to_treatment<=5 & treatment == 1 
	replace wm_outcome = 1 if arm == "NoTreatment" & time_to_treatment<=5 & treatment == 1 

	** Case 2: they deviate at 5 days
	replace wm_fup = 5 if arm == "NoTreatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 
	replace wm_outcome = 0 if arm == "NoTreatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "NoTreatment" & treatment == 0 & fup < 5 
	replace wm_outcome = 0 if arm == "NoTreatment" & treatment == 0 & fup < 5

	** add 1 days to 0-survivors
	replace wm_fup= 1 if arm == "NoTreatment" & wm_fup==0

	order pseudo_key arm
	compress
	save "matched_10/subgp`k'/cloned/cloned `event' subgp_`k' bs`j'.dta", replace
}
}
}
*
* Split times
qui foreach event in covid_hosp MMMI vag_bleed preg_ht preeclampsia eclampsia HELLP preterm_birth unspec_comp abx ICU HDU maternal_death c_section preterm_birth_comb stillbirth nnd {
forvalues k = 1/7 {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_`k'" _col(40) "bs`j'" _col(55) "split" _col(70) "matched_10"
	}
	use "matched_10/subgp`k'/cloned/cloned `event' subgp_`k' bs`j'.dta", clear

	** times of event
	stset fup, fail(outcome) id(NewID)
	stsplit, at(failures)

	gen tstart = _t0
	gen tstop = _t

	** times of censoring
	gen TrialEmul_cens = 1-outcome
	stset fup, fail(TrialEmul_cens) id(NewID)
	stsplit, at(failures)

	replace tstart = _t0 if tstart<_t0 & _t0 != . & _t != .
	replace tstop = _t if tstop>_t & _t0 != . & _t != .

	order pseudo_key arm tstart tstop
	sort _all
	compress

	save "matched_10/subgp`k'/split/split `event' subgp_`k' bs`j'.dta", replace
}
}
}
*
* IPCW - Treatment arm
qui foreach event in covid_hosp MMMI vag_bleed preg_ht preeclampsia eclampsia HELLP preterm_birth unspec_comp abx ICU HDU maternal_death c_section preterm_birth_comb stillbirth nnd {
forvalues k = 1/7 {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_`k'" _col(40) "bs`j'" _col(55) "Treatment" _col(70) "matched_10"
	}
	use "matched_10/subgp`k'/split/split `event' subgp_`k' bs`j'.dta", clear
	keep if arm == "Treatment"
	replace outcome = . if TrialEmul_cens == .
	merge m:1 pseudo_key using "covid pregnant antiviral outpatient characteristics.dta", keep(3) nogen

	* adapt records to the long format
	sort NewID tstart
	bysort NewID (tstart): replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	capture stcox c.age#i.age_gp i.vaccine_status i.month_covid i.immuno_hist i.HELLP_bl i.GDM i.MMMI_bl i.preg_ht_preeclampsia_bl, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "matched_10/subgp`k'/Treatment/Treatment `event' subgp_`k' bs`j'.dta", replace
}
}
}
*
* IPCW - Control arm
qui foreach event in covid_hosp MMMI vag_bleed preg_ht preeclampsia eclampsia HELLP preterm_birth unspec_comp abx ICU HDU maternal_death c_section preterm_birth_comb stillbirth nnd {
forvalues k = 1/7 {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_`k'" _col(40) "bs`j'" _col(55) "NoTreatment" _col(70) "matched_10"
	}
	use "matched_10/subgp`k'/split/split `event' subgp_`k' bs`j'.dta", clear
	keep if arm == "NoTreatment"
	replace outcome = . if TrialEmul_cens == .
	merge m:1 pseudo_key using "covid pregnant antiviral outpatient characteristics.dta", keep(3) nogen

	* adapt records to the long format
	sort NewID tstart
	bysort NewID: replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	capture stcox c.age#i.age_gp i.vaccine_status i.month_covid i.immuno_hist i.HELLP_bl i.GDM i.MMMI_bl i.preg_ht_preeclampsia_bl, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "matched_10/subgp`k'/NoTreatment/Notreatment `event' subgp_`k' bs`j'.dta", replace
}
}
}
*
* Combine & Generate weights
qui foreach event in covid_hosp MMMI vag_bleed preg_ht preeclampsia eclampsia HELLP preterm_birth unspec_comp abx ICU HDU maternal_death c_section preterm_birth_comb stillbirth nnd {
qui forvalues k = 1/7 {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_`k'" _col(40) "bs`j'" _col(55) "Combine" _col(70) "matched_10"
	}
	use "matched_10/subgp`k'/Treatment/Treatment `event' subgp_`k' bs`j'.dta", clear
	append using "matched_10/subgp`k'/NoTreatment/Notreatment `event' subgp_`k' bs`j'.dta"

	// create a new ID variable for each clone in each arm
	tostring NewID, gen(NewID_str)
	replace NewID_str = "00000" + NewID_str if length(NewID_str)== 1
	replace NewID_str = "0000" + NewID_str if length(NewID_str)== 2
	replace NewID_str = "000" + NewID_str if length(NewID_str)== 3
	replace NewID_str = "00" + NewID_str if length(NewID_str)== 4
	replace NewID_str = "0" + NewID_str if length(NewID_str)== 5
	gen Anal_ID = "1" + NewID_str if arm == "Treatment" 
	replace Anal_ID = "2" + NewID_str if arm == "NoTreatment" 

	replace weight = 1 if wm_outcome == 1 & tstop == 5

	gen arm_value = 1 if arm == "Treatment"
	replace arm_value = 0 if arm_value == .

	keep pseudo_key tstart tstop fup_obs event time_to_treatment treatment bs outcome fup wm_fup wm_outcome TrialEmul_cens weight Anal_ID arm_value invalid
	destring Anal_ID, replace
	compress
	save "matched_10/subgp`k'/emulated/emulated `event' subgp_`k' bs`j'.dta", replace
}
}
}
*
* Generate weights
qui foreach event in covid_hosp MMMI vag_bleed preg_ht preeclampsia eclampsia HELLP preterm_birth unspec_comp abx ICU HDU maternal_death c_section preterm_birth_comb stillbirth nnd {
forvalues k = 1/7 {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_`k'" _col(40) "bs`j'" _col(55) "ipcw" _col(70) "matched_10"
	}
	use "matched_10/subgp`k'/emulated/emulated `event' subgp_`k' bs`j'.dta", clear
	rename weight _ipcw
	* generate trimmed IPCW
	gen _ipcw_trim = _ipcw
	sum fup_obs, d
	forvalues t = 0/`r(max)' {
		forvalues i = 0/1 {
			sum _ipcw_trim if arm_value == `i' & tstart == `t', d
			replace _ipcw_trim = r(p1) if arm_value == `i' & tstart == `t' & _ipcw_trim < r(p1)
			replace _ipcw_trim = r(p99) if arm_value == `i' & tstart == `t' & _ipcw_trim > r(p99)
		}
	}
	*
	compress
	save "matched_10/subgp`k'/weight/weight `event' subgp_`k' bs`j'.dta", replace
}
}
}
*
cls
* Generate KM estimate
qui foreach event in covid_hosp MMMI vag_bleed preg_ht preeclampsia eclampsia HELLP preterm_birth unspec_comp abx ICU HDU maternal_death c_section preterm_birth_comb stillbirth nnd {
forvalues k = 1/7 {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_`k'" _col(40) "bs`j'" _col(55) "KM ipcw" _col(80) "matched_10"
	}
	use "matched_10/subgp`k'/weight/weight `event' subgp_`k' bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = _ipcw], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(fup bs invalid)
	save "matched_10/subgp`k'/KM ipcw/KM `event' subgp_`k' bs`j'.dta", replace
}
}
}
*
* Finalize bootstrap datasets
foreach event in covid_hosp MMMI vag_bleed preg_ht preeclampsia eclampsia HELLP preterm_birth unspec_comp abx ICU HDU maternal_death c_section preterm_birth_comb stillbirth nnd {
forvalues k = 1/7 {
	clear
forvalues j = 0/500 {
	append using "matched_10/subgp`k'/KM ipcw/KM `event' subgp_`k' bs`j'.dta"
}
	gen hazard_s_w = 1 - KM_s_w
	gen hazard_ns_w = 1 - KM_ns_w
	gen odds_s_w = hazard_s_w/(1-hazard_s_w)
	gen odds_ns_w = hazard_ns_w/(1-hazard_ns_w)
	gen RR_w = hazard_s_w/hazard_ns_w
	gen diff_w = hazard_s_w - hazard_ns_w
	gen OR_w = odds_s_w / odds_ns_w
	compress
	save "matched_10/subgp`k'/KM ipcw `event' subgp_`k' bs_all.dta", replace
}
}
*

**# Table 2 - Main results
cls
* N of events
qui forvalues j = 1/1 {
	if `j' == 1 {
	noi di "event" _col(15) "subgp" _col(30) "N_1" _col(45) "N_0"
}
foreach event in covid_hosp MMMI vag_bleed preg_ht preeclampsia eclampsia HELLP preterm_birth unspec_comp abx ICU HDU maternal_death c_section preterm_birth_comb stillbirth nnd {
forvalues k = 1/1 {
	use "matched_10/subgp`k'/prepare/prepare `event' subgp_`k' bs0.dta", clear
	count if event == 1 & treatment == 1
	scalar N_1 = r(N)
	count if event == 1 & treatment == 0
	scalar N_0 = r(N)
	noi di substr("`event'",1,13) _col(15) "subgp_`k'" _col(30) N_1 _col(45) N_0
}
}
}
*
* Cumulative incidence (N, %)
qui foreach event in covid_hosp MMMI vag_bleed preg_ht preeclampsia eclampsia HELLP preterm_birth unspec_comp abx ICU HDU maternal_death c_section preterm_birth_comb stillbirth nnd {
forvalues k = 1/1 {
	use "matched_10/subgp`k'/KM ipcw `event' subgp_`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum hazard_s_w if bs == 0, d
	scalar hazard_s_mean = r(mean)
	centile hazard_s_w if bs > 0, centile(2.5 97.5)
	scalar hazard_s_cil = r(c_1)
	scalar hazard_s_ciu = r(c_2)
	
	sum hazard_ns_w if bs == 0, d
	scalar hazard_ns_mean = r(mean)
	centile hazard_ns_w if bs > 0, centile(2.5 97.5)
	scalar hazard_ns_cil = r(c_1)
	scalar hazard_ns_ciu = r(c_2)
	
	noi di substr("`event'",1,13) _col(15) "subgp_`k'" _col(30) hazard_s_mean _col(45) hazard_ns_mean
}
}
*
* Absolute risk reduction
qui foreach event in covid_hosp MMMI vag_bleed preg_ht preeclampsia eclampsia HELLP preterm_birth unspec_comp abx ICU HDU maternal_death c_section preterm_birth_comb stillbirth nnd {
forvalues k = 1/1 {
	use "matched_10/subgp`k'/KM ipcw `event' subgp_`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum diff_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile diff_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "subgp_`k'" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
}
*
* Relative risk
qui foreach event in covid_hosp MMMI vag_bleed preg_ht preeclampsia eclampsia HELLP preterm_birth unspec_comp abx ICU HDU maternal_death c_section preterm_birth_comb stillbirth nnd {
qui forvalues k = 1/1 {
	use "matched_10/subgp`k'/KM ipcw `event' subgp_`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum RR_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	sum RR_w, d
	centile RR_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "subgp_`k'" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
}
*

********************************************************************************
********************************************************************************
**# Supplementary Figure 1 - PS distribution plot
* before 1:10 PS matching
capture mkdir "graph"
use "covid pregnant antiviral outpatient matched subgp_1.dta", clear
bysort group_old (prob_group) : replace prob_group = int(prob_group/0.02)*0.02 if _n == 1
twoway (hist prob_group if group_old == 1, lcolor(navy) fcolor(navy%40) width(0.02) percent) ///
(hist prob_group if group_old == 0, lpattern(dash) lcolor(maroon) fcolor(maroon%40) width(0.02) percent), ///
title("Propensity score distribution before matching") ///
xtitle("Propensity score") ytitle("Proportion (%)") ///
legend( label( 1 "Nirmatrelvir/ritonavir") label(2 "Control") ) ///
plotregion(color(white)) graphregion(color(white))  ///
ylabel(, nogrid) xlabel(0(0.1)1, nogrid)
graph export "graph\ps_group_10_before.png", replace

* after 1:10 PS matching
use "covid pregnant antiviral outpatient matched subgp_1.dta", clear
bysort group_10 (prob_group) : replace prob_group = int(prob_group/0.02)*0.02 if _n == 1
twoway (hist prob_group if group_10 == 1, lcolor(navy) fcolor(navy%40) width(0.02) percent) ///
(hist prob_group if group_10 == 0, lpattern(dash) lcolor(maroon) fcolor(maroon%40) width(0.02) percent), ///
title("Propensity score distribution after matching") ///
xtitle("Propensity score") ytitle("Proportion (%)") ///
legend( label( 1 "Nirmatrelvir/ritonavir") label(2 "Control" ) ) ///
plotregion(color(white)) graphregion(color(white))  ///
ylabel(, nogrid) xlabel(0(0.1)1, nogrid)
graph export "graph\ps_group_10_after.png", replace

********************************************************************************
**# Supplementary Figure 2 - Distribution of IPCW
use "matched_10/subgp1/weight/weight covid_hosp subgp_1 bs0.dta", clear
keep if inrange(tstart,0,4)
replace tstart = tstart + 1
gen arm_value_0 = 1 - arm_value
gen ipcw_log = log10(_ipcw)
keep *_log arm_value_0 tstart

graph box ipcw_log, over(tstart, relabel(1 "Day 1" 2 "Day 2" 3 "Day 3" 4 "Day 4" 5 "Day 5")) over(arm_value_0, relabel(1 "Nirmatrelvir/ritonavir" 2 "Control")) ytitle("Logarithm of the IPCW")
graph export "graph\ipcw_log.png", replace

********************************************************************************
**# Supplementary Figure 3 - Cumulative incidence plot (28-day MMMI)
* MMMI
cls
use "matched_10/subgp1/weight/weight MMMI subgp_1 bs0.dta", clear
stset tstop, enter(time tstart) failure(outcome)
sts graph, failure by(arm_value) xlabel(0(7)28) ytitle("%", orientation(horizontal)) ///
ylabel(0 "0" .005 "0.5" .01 "1" .015 "1.5" .02 "2", format(%5.1f)) tmax(28) ttitle("Days") ///
plot1opts(lcolor(maroon)lpattern(dash)) plot2(lcolor(navy)) ///
plotregion(color(white)) graphregion(color(white)) ///
title("MMMI") ///
legend(order(2 "Nirmatrelvir/ritonavir" 1 "Control")) ///
risktable(,justification(center) size(*0.8) order(2 "Nirmatrelvir/ritonavir" 1 "Control") rowtitle(, justification(left)))
graph export "graph\KM_MMMI.png", replace

********************************************************************************
********************************************************************************
**# Supplementary Table 1 - Baseline characteristics before 1:10 ps matching
cls
* N/mean, %/SD
qui forvalues k = 1(-1)0 {
	use "covid pregnant antiviral outpatient matched subgp_1.dta", clear
	replace name_dose_1_num = . if name_dose_1_num == 4
	gen covid_infect0 = 1 - covid_infect
	keep if group_old == `k'
	noi di _newline "before"
	noisily di "group=`k'"
	noisily di "var" _col(25) "N" _col(35) "mean" _col(50) "sd"
	foreach var in age {
		sum `var', d
		scalar m1 = r(mean)
		scalar n1 = r(N)
		scalar sd1 = r(sd)
		noi di "`var'" _col(25) n1 _col(35) m1 _col(50) sd1
	}

	noisily di "var" _col(25) "N" _col(35) "%"
	foreach var in covid_infect0 covid_infect vaccine_full chd ckd liver cancer GDM cholestasis MMMI_bl vag_bleed_bl preg_ht_bl preeclampsia_bl eclampsia_bl HELLP_bl unspec_comp_bl {
		sum `var', d
		scalar m1 = r(mean)
		scalar n1 = r(N)
		scalar b = round(n1*m1)
		noisily di "`var'" _col(25) b _col(35) m1*100
	}
	foreach var in age_gp vaccine_status_cate {
		noi tab `var'
	}
}
*
* SMD
qui forvalues j = 1/1 {
	use "covid pregnant antiviral outpatient matched subgp_`j'.dta", clear
	replace name_dose_1_num = . if name_dose_1_num == 4
	replace age_gp = age > 32
	noi di _newline "var" _col(25) "SMD"
***continuous
	foreach var in age {
		stddiff `var', by(group_old) abs cohensd
		noisily di "`var'" _col(25) abs(r(stddiff)[1,1])
	}
	noi di
***categorical
* binary
	foreach var in covid_infect vaccine_full ///
	chd ckd liver cancer GDM cholestasis MMMI_bl vag_bleed_bl preg_ht_bl preeclampsia_bl eclampsia_bl HELLP_bl unspec_comp_bl ///
	age_gp vaccine_status_cate {
		tab `var' group_old
		if r(r) > 1 { 
			stddiff i.`var', by(group_old) abs cohensd
			noisily di "`var'" _col(25) abs(r(stddiff)[1,1])
		}
		else {
			noisily di "`var'" _col(25) .
		}
	}
	noisily di
}
*

********************************************************************************
**# Supplementary Table 2 - Baseline characteristics (5-day grace period) after 1:10 ps matching and IPCW
forvalues t = 1/5 {
	use "matched_10/subgp1/weight/weight covid_hosp subgp_1 bs0.dta", clear
	merge m:1 pseudo_key using "covid pregnant antiviral outpatient matched subgp_1", keep(3) nogen
	keep if tstart == `t' - 1
	compress
	save "matched_10/emulated covid_hosp subgp_1 bs0 day`t'", replace
}
*
cls
* N/mean, %/SD
qui forvalues t = 1/5 {
forvalues k = 1(-1)0 {
	use "matched_10/emulated covid_hosp subgp_1 bs0 day`t'", clear
	replace arm_value = . if date_death <= date_baseline
	gen gender0 = 1 - gender
	gen number_dose_cate = number_dose >= 3
	gen number_dose_cate_0 = number_dose_cate == 0
	replace name_dose_1_num = . if name_dose_1_num == 4
	gen covid_infect0 = 1 - covid_infect
	
	noi di _newline "after"
	noisily di "arm_value=`k', day`t'"
	
	***continuous variables
	noisily di "var" _col(30) "N" _col(45) "mean" _col(60) "sd"
	foreach var in age charlson_index onset_dur gest_wk {
		sum `var' [w=_ipcw] if arm_value == `k'
		noi di "`var'" _col(30) r(N) _col(45) r(mean) _col(60) r(sd)
	}
	*** categorical variables
	noisily di "var" _col(30) "N" _col(45) "%"
	foreach var in gender gender0 covid_infect0 covid_infect immuno_hist vaccine_full number_dose_cate_0 number_dose_cate dm ht lung chd ckd liver cancer gout obesity obs_sleep_apnoea depression GDM cholestasis MMMI_bl vag_bleed_bl preg_ht_bl preeclampsia_bl preg_ht_preeclampsia_bl eclampsia_bl HELLP_bl unspec_comp_bl hypothyroidism hypoglycaemia drug_acei_arb drug_anticoagulant drug_antiplatelet drug_lipid drug_nsaid drug_bb drug_ccb drug_diuretics drug_antidepressants {
		sum `var' if arm_value == `k' [w=_ipcw]
		noi di "`var'=1" _col(30) r(N)*r(mean) _col(45) r(mean)*100
	}
	foreach var in age_gp cci_gp vaccine_status_cate gest_wk_cate {
	forvalues j = 1/2 {
		gen `var'`j' = `var' == `j'
		sum `var'`j' if arm_value == `k' [w=_ipcw]
		noi di "`var'=`j'" _col(30) r(N)*r(mean) _col(45) r(mean)*100
	}
}
}
}
*
* SMD
qui forvalues t = 1/5 {
	noi di "day"`t'
forvalues j = 1/1 {
	use "matched_10/emulated covid_hosp subgp_1 bs0 day`t'", clear
	replace gest_wk = 0 if gest_wk == .
	replace arm_value = . if date_death <= date_baseline
	gen number_dose_cate = number_dose >= 3
	gen number_dose_cate_0 = number_dose_cate == 0
	replace name_dose_1_num = . if name_dose_1_num == 4
	replace age_gp = age > 32
	noi di "var" _col(30) "SMD"
	
	***continuous variables
	foreach var in age charlson_index onset_dur gest_wk {
		sum `var' if arm_value == 1 [w=_ipcw]
		local mean_1 = r(mean)
		local sd_1 = r(sd)
		sum `var' if arm_value == 0 [w=_ipcw]
		local mean_0 = r(mean)
		local sd_0 = r(sd)
		stddiffi `mean_1' `sd_1' `mean_0' `sd_0'
		noi di "`var'" _col(30) abs(r(std_diff))
	}
	noi di
	*** categorical variables
	foreach var in gender covid_infect immuno_hist vaccine_full number_dose_cate dm ht lung chd ckd liver cancer gout obesity obs_sleep_apnoea depression GDM cholestasis MMMI_bl vag_bleed_bl preg_ht_bl preeclampsia_bl preg_ht_preeclampsia_bl eclampsia_bl HELLP_bl unspec_comp_bl hypothyroidism hypoglycaemia drug_acei_arb drug_anticoagulant drug_antiplatelet drug_lipid drug_nsaid drug_bb drug_ccb drug_diuretics drug_antidepressants age_gp cci_gp vaccine_status_cate gest_wk_cate {
	tab `var' arm_value
	if r(r) == 2 {
		tab `var' if arm_value == 1 [aw=_ipcw], matcell(temp)
		local m11 = int(temp[1,1])
		local m21 = int(temp[2,1])
		tab `var' if arm_value == 0 [aw=_ipcw], matcell(temp)
		local m12 = int(temp[1,1])
		local m22 = int(temp[2,1])
		capture stddiffi `m11' `m12' \ `m21' `m22'
	}
	if r(r) == 3 {
		tab `var' if arm_value == 1 [aw=_ipcw], matcell(temp)
		local m11 = int(temp[1,1])
		local m21 = int(temp[2,1])
		local m31 = int(temp[3,1])
		tab `var' if arm_value == 0 [aw=_ipcw], matcell(temp)
		local m12 = int(temp[1,1])
		local m22 = int(temp[2,1])
		local m32 = int(temp[3,1])
		capture stddiffi `m11' `m12' \ `m21' `m22' \ `m31' `m32'
	}
	noi di "`var'" _col(30) abs(r(std_diff))
}
}
}
*

********************************************************************************
**# Supplementary Table 3 - Secondary diagnosis of hospitalization among nirmatrelvir/ritonavir and control groups
use "dx_extracted.dta", clear
merge m:1 pseudo_key using "date_covid_hosp.dta", keep(3) keepusing(pseudo_key) nogen
merge m:1 pseudo_key using "covid pregnant antiviral outpatient matched subgp_1.dta", keep(3) keepusing(date_baseline group_10) nogen
keep if reference_date_id > date_baseline
keep if group_10 < .

tab icd9_cd group_10

* drop primary reason (SARS-CoV-2 infection)
drop if inlist(icd9_cd, "079.89", "519.8")

tab icd9_cd group_10

********************************************************************************
**# Supplementary Table 4 - Subgroups analysis
cls
* Cumulative incidence (N, %)
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues k = 1/7 {
	use "matched_10/subgp`k'/KM ipcw `event' subgp_`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum hazard_s_w if bs == 0, d
	scalar hazard_s_mean = r(mean)
	centile hazard_s_w if bs > 0, centile(2.5 97.5)
	scalar hazard_s_cil = r(c_1)
	scalar hazard_s_ciu = r(c_2)
	
	sum hazard_ns_w if bs == 0, d
	scalar hazard_ns_mean = r(mean)
	centile hazard_ns_w if bs > 0, centile(2.5 97.5)
	scalar hazard_ns_cil = r(c_1)
	scalar hazard_ns_ciu = r(c_2)
	
	noi di substr("`event'",1,13) _col(15) "subgp_`k'" _col(30) hazard_s_mean _col(45) hazard_ns_mean
}
}
*
* Absolute risk reduction
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues k = 1/7 {
	use "matched_10/subgp`k'/KM ipcw `event' subgp_`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum diff_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile diff_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "subgp_`k'" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
}
*
* Relative risk
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
qui forvalues k = 1/7 {
	use "matched_10/subgp`k'/KM ipcw `event' subgp_`k' bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum RR_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	sum RR_w, d
	centile RR_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "subgp_`k'" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
}
*

********************************************************************************
**# Supplementary Table 5 - Sensitivity analyses

capture mkdir "Sensitivity"

*** Sen1 - Excluding patients whose index date was defined by nirmatrelvir/ritonavir initiation

***Creates directory
capture mkdir "Sensitivity/SEN1"
capture mkdir "Sensitivity/SEN1/prepare"
capture mkdir "Sensitivity/SEN1/cloned"
capture mkdir "Sensitivity/SEN1/split"
capture mkdir "Sensitivity/SEN1/Treatment"
capture mkdir "Sensitivity/SEN1/NoTreatment"
capture mkdir "Sensitivity/SEN1/emulated"
capture mkdir "Sensitivity/SEN1/weight"
capture mkdir "Sensitivity/SEN1/KM ipcw"
*

*** Define SEN1
use "covid pregnant antiviral outpatient.dta", clear
gen date_baseline_report = min(date_covid_report, date_onset)
drop if date_baseline_report != date_baseline
drop date_baseline_report

gen covid_infect = (date_covid >= mdy(9,1,2022))

local mi = "logit group age c.date_covid#i.covid_infect i.vaccine_status age##i.MMMI_bl i.immuno_hist i.HELLP_bl i.GDM i.MMMI_bl i.preg_ht_preeclampsia_bl"
	
`mi'
predict prob_group
gen xb_group = ln(prob_group / (1-prob_group))

rename group group_old
noi di "matching 1:10"
set seed 123456
calipmatch, gen(group_10_match) case(group_old) max(10) calipermatch(xb_group) caliperwidth(.05)
bysort pseudo_key (group_10_match): replace group_10_match = group_10_match[1]
gen group_10 = group_old if group_10_match < .

compress
save "Sensitivity/SEN1/covid pregnant antiviral outpatient matched SEN1.dta", replace

* Prepare dataset for analysis
cls
qui foreach event in covid_hosp MMMI {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "setup" _col(70) "SEN1"
	}
	* Setup for trial emulation
	use "Sensitivity/SEN1/covid pregnant antiviral outpatient matched SEN1.dta", clear
	rename group_10 group
	keep if group < .
	
	* events
	merge 1:1 pseudo_key using "date_covid_hosp", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_MMMI", keep(1 3) nogen
	
	if `j' > 0 {
		set seed `j'
		bsample
	}
	replace date_paxlovid = . if date_paxlovid >= date_admission
	* organize / rename / generate variables
	gen date_last_fu = min(date_death, mdy(02,12,2023), date_baseline + 28)
	gen date_event = date_`event'
	gen event = inrange(date_event, date_baseline, date_last_fu)
	gen fup_obs = min(date_event-date_baseline, 28) if event == 1
	replace fup_obs = min(date_last_fu-date_baseline, 28) if event == 0
	tab fup_obs
	gen time_to_treatment = date_paxlovid - date_baseline
	tab time_to_treatment
	gen treatment = inrange(date_paxlovid - date_baseline, 0, 5)
	* keep necessary variables
	keep pseudo_key fup_obs event time_to_treatment treatment
	gen bs = `j'
	compress
	save "Sensitivity/SEN1/prepare/prepare `event' SEN1 bs`j'.dta", replace
}
}
*
qui foreach event in c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "setup" _col(70) "SEN1"
	}
	* Setup for trial emulation
	use "Sensitivity/SEN1/covid pregnant antiviral outpatient matched SEN1.dta", clear
	rename group_10 group
	keep if group < .
	
	* events
	merge 1:1 pseudo_key using "date_c_section", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_preterm_birth", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_stillbirth", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_neonatal_death", keep(1 3) nogen
	
	if `j' > 0 {
		set seed `j'
		bsample
	}
	replace date_paxlovid = . if date_paxlovid >= date_admission
	* organize / rename / generate variables
	gen date_last_fu = min(date_death, mdy(02,12,2023))
	gen date_event = date_`event'
	gen event = inrange(date_event, date_baseline, date_last_fu)
	gen fup_obs = date_event-date_baseline if event == 1
	replace fup_obs = date_last_fu-date_baseline if event == 0
	tab fup_obs
	gen time_to_treatment = date_paxlovid - date_baseline
	tab time_to_treatment
	gen treatment = inrange(date_paxlovid - date_baseline, 0, 5)
	* keep necessary variables
	keep pseudo_key fup_obs event time_to_treatment treatment
	gen bs = `j'
	compress
	save "Sensitivity/SEN1/prepare/prepare `event' SEN1 bs`j'.dta", replace
}
}
*
cls
***** Bootstrap
*** Cloning & censoring
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "cloning" _col(70) "SEN1" 
	}
	* Prepare dataset for analysis
	use "Sensitivity/SEN1/prepare/prepare `event' SEN1 bs`j'.dta", replace
	stset fup_obs, failure(event)

	* Arm A: no treatment within 5 days (control: non-exposed group)
	gen outcomeA = _d // _d = `event'
	gen fupA = _t // _t = follow up time

	/// if the patient received treatment within 5 days:
	/// 1. no event outcome, since the patient survived till censoring (treatment)
	replace outcomeA = 0 if treatment==1 & time_to_treatment <=5 
	/// 2. follow up is censored at treatment
	replace fupA = time_to_treatment if treatment==1 & time_to_treatment <=5

	* Arm B: treatment within 5 days (treated: exposed group)
	gen outcomeB = _d 
	gen fupB = _t 

	/// if the patient survived the first 5 days and did not receive treatment within 5 days:
	/// 1. no event outcome if the patient survived the first 5 days
	replace outcomeB = 0 if (treatment==0 & _t>5) | (treatment==1 & time_to_treatment >5 & time_to_treatment !=.)
	/// 2. follow up is censored at 5 days
	replace fupB = 5 if (treatment==0 & _t>5) | (treatment==1 & time_to_treatment >5 & time_to_treatment != .)

	** append clones 
	preserve
		drop outcomeB fupB
		rename outcomeA outcome
		rename fupA fup
		gen arm = "NoTreatment"
		tempfile a
		save "`a'", replace
	restore
		drop outcomeA fupA
		rename outcomeB outcome
		rename fupB fup
		gen arm = "Treatment"	
		cap append using "`a'"

	// Weight models

	sort _all
	gen NewID = _n

	** add 1 day to 0-survivors
	replace fup= 1 if fup==0

	** Weight model: define survival time and event indicator	
	* treatment Arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	gen wm_fup = time_to_treatment if arm == "Treatment" & time_to_treatment<=5 & time_to_treatment!=. & treatment == 1 
	gen wm_outcome = 0 if arm == "Treatment" & time_to_treatment<=5 & time_to_treatment!=. & treatment == 1 

	** Case 2: they deviate at 5 days
	replace wm_fup = 5 if arm == "Treatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1))
	replace wm_outcome = 1 if arm == "Treatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "Treatment" & treatment == 0 & fup < 5
	replace wm_outcome = 0 if arm == "Treatment" & treatment == 0 & fup < 5
	** add 1 days to 0-survivors
	replace wm_fup= 1 if arm == "Treatment" & wm_fup==0 

	* No treatment Arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	replace wm_fup = time_to_treatment if arm == "NoTreatment" & time_to_treatment<=5 & treatment == 1 
	replace wm_outcome = 1 if arm == "NoTreatment" & time_to_treatment<=5 & treatment == 1 

	** Case 2: they deviate at 5 days
	replace wm_fup = 5 if arm == "NoTreatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 
	replace wm_outcome = 0 if arm == "NoTreatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "NoTreatment" & treatment == 0 & fup < 5 
	replace wm_outcome = 0 if arm == "NoTreatment" & treatment == 0 & fup < 5

	** add 1 days to 0-survivors
	replace wm_fup= 1 if arm == "NoTreatment" & wm_fup==0

	order pseudo_key arm
	compress
	save "Sensitivity/SEN1/cloned/cloned `event' SEN1 bs`j'.dta", replace
}
}
*
* Split times
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "split" _col(70) "SEN1"
	}
	use "Sensitivity/SEN1/cloned/cloned `event' SEN1 bs`j'.dta", clear

	** times of event
	stset fup, fail(outcome) id(NewID)
	stsplit, at(failures)

	gen tstart = _t0
	gen tstop = _t

	** times of censoring
	gen TrialEmul_cens = 1-outcome
	stset fup, fail(TrialEmul_cens) id(NewID)
	stsplit, at(failures)

	replace tstart = _t0 if tstart<_t0 & _t0 != . & _t != .
	replace tstop = _t if tstop>_t & _t0 != . & _t != .

	order pseudo_key arm tstart tstop
	sort _all
	compress

	save "Sensitivity/SEN1/split/split `event' SEN1 bs`j'.dta", replace
}
}
*
* IPCW - Treatment arm
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "Treatment" _col(70) "SEN1"
	}
	use "Sensitivity/SEN1/split/split `event' SEN1 bs`j'.dta", clear
	keep if arm == "Treatment"
	replace outcome = . if TrialEmul_cens == .
	merge m:1 pseudo_key using "covid pregnant antiviral outpatient characteristics.dta", keep(3) nogen

	* adapt records to the long format
	sort NewID tstart
	bysort NewID (tstart): replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	capture stcox c.age#i.age_gp i.vaccine_status i.month_covid i.immuno_hist i.HELLP_bl i.GDM i.MMMI_bl i.preg_ht_preeclampsia_bl, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "Sensitivity/SEN1/Treatment/Treatment `event' SEN1 bs`j'.dta", replace
}
}
*
* IPCW - Control arm
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "NoTreatment" _col(70) "SEN1"
	}
	use "Sensitivity/SEN1/split/split `event' SEN1 bs`j'.dta", clear
	keep if arm == "NoTreatment"
	replace outcome = . if TrialEmul_cens == .
	merge m:1 pseudo_key using "covid pregnant antiviral outpatient characteristics.dta", keep(3) nogen

	* adapt records to the long format
	sort NewID tstart
	bysort NewID: replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	capture stcox c.age#i.age_gp i.vaccine_status i.month_covid i.immuno_hist i.HELLP_bl i.GDM i.MMMI_bl i.preg_ht_preeclampsia_bl, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "Sensitivity/SEN1/NoTreatment/Notreatment `event' SEN1 bs`j'.dta", replace
}
}
*
* Combine & Generate weights
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "Combine" _col(70) "SEN1"
	}
	use "Sensitivity/SEN1/Treatment/Treatment `event' SEN1 bs`j'.dta", clear
	append using "Sensitivity/SEN1/NoTreatment/Notreatment `event' SEN1 bs`j'.dta"

	// create a new ID variable for each clone in each arm
	tostring NewID, gen(NewID_str)
	replace NewID_str = "00000" + NewID_str if length(NewID_str)== 1
	replace NewID_str = "0000" + NewID_str if length(NewID_str)== 2
	replace NewID_str = "000" + NewID_str if length(NewID_str)== 3
	replace NewID_str = "00" + NewID_str if length(NewID_str)== 4
	replace NewID_str = "0" + NewID_str if length(NewID_str)== 5
	gen Anal_ID = "1" + NewID_str if arm == "Treatment" 
	replace Anal_ID = "2" + NewID_str if arm == "NoTreatment" 

	replace weight = 1 if wm_outcome == 1 & tstop == 5

	gen arm_value = 1 if arm == "Treatment"
	replace arm_value = 0 if arm_value == .

	keep pseudo_key tstart tstop fup_obs event time_to_treatment treatment bs outcome fup wm_fup wm_outcome TrialEmul_cens weight Anal_ID arm_value invalid
	destring Anal_ID, replace
	compress
	save "Sensitivity/SEN1/emulated/emulated `event' SEN1 bs`j'.dta", replace
}
}
*
* Generate weights
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "ipcw" _col(70) "SEN1"
	}
	use "Sensitivity/SEN1/emulated/emulated `event' SEN1 bs`j'.dta", clear
	rename weight _ipcw
	compress
	save "Sensitivity/SEN1/weight/weight `event' SEN1 bs`j'.dta", replace
}
}
*
cls
* Generate KM estimate
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "KM ipcw" _col(80) "SEN1"
	}
	use "Sensitivity/SEN1/weight/weight `event' SEN1 bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = _ipcw], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(fup bs invalid)
	save "Sensitivity/SEN1/KM ipcw/KM `event' SEN1 bs`j'.dta", replace
}
}
*
* Finalize bootstrap datasets
foreach event in covid_hosp MMMI c_section preterm_birth_comb {
clear
forvalues j = 0/500 {
	append using "Sensitivity/SEN1/KM ipcw/KM `event' SEN1 bs`j'.dta"
}
	gen hazard_s_w = 1 - KM_s_w
	gen hazard_ns_w = 1 - KM_ns_w
	gen odds_s_w = hazard_s_w/(1-hazard_s_w)
	gen odds_ns_w = hazard_ns_w/(1-hazard_ns_w)
	gen RR_w = hazard_s_w/hazard_ns_w
	gen diff_w = hazard_s_w - hazard_ns_w
	gen OR_w = odds_s_w / odds_ns_w
	compress
	save "Sensitivity/SEN1/KM ipcw `event' SEN1 bs_all.dta", replace
}
*
cls
* Cumulative incidence (N, %)
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	use "Sensitivity/SEN1/KM ipcw `event' SEN1 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum hazard_s_w if bs == 0, d
	scalar hazard_s_mean = r(mean)
	centile hazard_s_w if bs > 0, centile(2.5 97.5)
	scalar hazard_s_cil = r(c_1)
	scalar hazard_s_ciu = r(c_2)
	
	sum hazard_ns_w if bs == 0, d
	scalar hazard_ns_mean = r(mean)
	centile hazard_ns_w if bs > 0, centile(2.5 97.5)
	scalar hazard_ns_cil = r(c_1)
	scalar hazard_ns_ciu = r(c_2)
	
	noi di substr("`event'",1,13) _col(15) "SEN1" _col(30) hazard_s_mean _col(45) hazard_ns_mean
}
*
* Absolute risk reduction
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	use "Sensitivity/SEN1/KM ipcw `event' SEN1 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum diff_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile diff_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "SEN1" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
*
* Relative risk
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	use "Sensitivity/SEN1/KM ipcw `event' SEN1 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum RR_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)

	centile RR_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "SEN1" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
*

******************************************
*** Sen2 - Excluding patients who were initiated nirmatrelvir/ritonavir beyond five days from the index date

***Creates directory
capture mkdir "Sensitivity/SEN2"
capture mkdir "Sensitivity/SEN2/prepare"
capture mkdir "Sensitivity/SEN2/cloned"
capture mkdir "Sensitivity/SEN2/split"
capture mkdir "Sensitivity/SEN2/Treatment"
capture mkdir "Sensitivity/SEN2/NoTreatment"
capture mkdir "Sensitivity/SEN2/emulated"
capture mkdir "Sensitivity/SEN2/weight"
capture mkdir "Sensitivity/SEN2/KM ipcw"
*

*** Define SEN2
use "covid pregnant antiviral outpatient.dta", clear
drop if date_paxlovid - date_baseline > 5 & date_paxlovid < .

gen covid_infect = (date_covid >= mdy(9,1,2022))

local mi = "logit group age c.date_covid#i.covid_infect i.vaccine_status age##i.MMMI_bl i.immuno_hist i.HELLP_bl i.GDM i.MMMI_bl i.preg_ht_preeclampsia_bl"
	
`mi'
predict prob_group
gen xb_group = ln(prob_group / (1-prob_group))

rename group group_old
noi di "matching 1:10"
set seed 123456
calipmatch, gen(group_10_match) case(group_old) max(10) calipermatch(xb_group) caliperwidth(.05)
bysort pseudo_key (group_10_match): replace group_10_match = group_10_match[1]
gen group_10 = group_old if group_10_match < .

compress
save "Sensitivity/SEN2/covid pregnant antiviral outpatient matched SEN2.dta", replace

* Prepare dataset for analysis
cls
qui foreach event in covid_hosp MMMI {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "setup" _col(70) "SEN2"
	}
	* Setup for trial emulation
	use "Sensitivity/SEN2/covid pregnant antiviral outpatient matched SEN2.dta", clear
	rename group_10 group
	keep if group < .
	
	* events
	merge 1:1 pseudo_key using "date_covid_hosp", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_MMMI", keep(1 3) nogen
	
	if `j' > 0 {
		set seed `j'
		bsample
	}
	replace date_paxlovid = . if date_paxlovid >= date_admission
	* organize / rename / generate variables
	gen date_last_fu = min(date_death, mdy(02,12,2023), date_baseline + 28)
	gen date_event = date_`event'
	gen event = inrange(date_event, date_baseline, date_last_fu)
	gen fup_obs = min(date_event-date_baseline, 28) if event == 1
	replace fup_obs = min(date_last_fu-date_baseline, 28) if event == 0
	tab fup_obs
	gen time_to_treatment = date_paxlovid - date_baseline
	tab time_to_treatment
	gen treatment = inrange(date_paxlovid - date_baseline, 0, 5)
	* keep necessary variables
	keep pseudo_key fup_obs event time_to_treatment treatment
	gen bs = `j'
	compress
	save "Sensitivity/SEN2/prepare/prepare `event' SEN2 bs`j'.dta", replace
}
}
*
qui foreach event in c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "setup" _col(70) "SEN2"
	}
	* Setup for trial emulation
	use "Sensitivity/SEN2/covid pregnant antiviral outpatient matched SEN2.dta", clear
	rename group_10 group
	keep if group < .
	
	* events
	merge 1:1 pseudo_key using "date_c_section", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_preterm_birth", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_stillbirth", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_neonatal_death", keep(1 3) nogen
	
	if `j' > 0 {
		set seed `j'
		bsample
	}
	replace date_paxlovid = . if date_paxlovid >= date_admission
	* organize / rename / generate variables
	gen date_last_fu = min(date_death, mdy(02,12,2023))
	gen date_event = date_`event'
	gen event = inrange(date_event, date_baseline, date_last_fu)
	gen fup_obs = date_event-date_baseline if event == 1
	replace fup_obs = date_last_fu-date_baseline if event == 0
	tab fup_obs
	gen time_to_treatment = date_paxlovid - date_baseline
	tab time_to_treatment
	gen treatment = inrange(date_paxlovid - date_baseline, 0, 5)
	* keep necessary variables
	keep pseudo_key fup_obs event time_to_treatment treatment
	gen bs = `j'
	compress
	save "Sensitivity/SEN2/prepare/prepare `event' SEN2 bs`j'.dta", replace
}
}
*
cls
***** Bootstrap
*** Cloning & censoring
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "cloning" _col(70) "SEN2" 
	}
	* Prepare dataset for analysis
	use "Sensitivity/SEN2/prepare/prepare `event' SEN2 bs`j'.dta", replace
	stset fup_obs, failure(event)

	* Arm A: no treatment within 5 days (control: non-exposed group)
	gen outcomeA = _d // _d = `event'
	gen fupA = _t // _t = follow up time

	/// if the patient received treatment within 5 days:
	/// 1. no event outcome, since the patient survived till censoring (treatment)
	replace outcomeA = 0 if treatment==1 & time_to_treatment <=5 
	/// 2. follow up is censored at treatment
	replace fupA = time_to_treatment if treatment==1 & time_to_treatment <=5

	* Arm B: treatment within 5 days (treated: exposed group)
	gen outcomeB = _d 
	gen fupB = _t 

	/// if the patient survived the first 5 days and did not receive treatment within 5 days:
	/// 1. no event outcome if the patient survived the first 5 days
	replace outcomeB = 0 if (treatment==0 & _t>5) | (treatment==1 & time_to_treatment >5 & time_to_treatment !=.)
	/// 2. follow up is censored at 5 days
	replace fupB = 5 if (treatment==0 & _t>5) | (treatment==1 & time_to_treatment >5 & time_to_treatment != .)

	** append clones 
	preserve
		drop outcomeB fupB
		rename outcomeA outcome
		rename fupA fup
		gen arm = "NoTreatment"
		tempfile a
		save "`a'", replace
	restore
		drop outcomeA fupA
		rename outcomeB outcome
		rename fupB fup
		gen arm = "Treatment"	
		cap append using "`a'"

	// Weight models

	sort _all
	gen NewID = _n

	** add 1 day to 0-survivors
	replace fup= 1 if fup==0

	** Weight model: define survival time and event indicator	
	* treatment Arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	gen wm_fup = time_to_treatment if arm == "Treatment" & time_to_treatment<=5 & time_to_treatment!=. & treatment == 1 
	gen wm_outcome = 0 if arm == "Treatment" & time_to_treatment<=5 & time_to_treatment!=. & treatment == 1 

	** Case 2: they deviate at 5 days
	replace wm_fup = 5 if arm == "Treatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1))
	replace wm_outcome = 1 if arm == "Treatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "Treatment" & treatment == 0 & fup < 5
	replace wm_outcome = 0 if arm == "Treatment" & treatment == 0 & fup < 5
	** add 1 days to 0-survivors
	replace wm_fup= 1 if arm == "Treatment" & wm_fup==0 

	* No treatment Arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	replace wm_fup = time_to_treatment if arm == "NoTreatment" & time_to_treatment<=5 & treatment == 1 
	replace wm_outcome = 1 if arm == "NoTreatment" & time_to_treatment<=5 & treatment == 1 

	** Case 2: they deviate at 5 days
	replace wm_fup = 5 if arm == "NoTreatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 
	replace wm_outcome = 0 if arm == "NoTreatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "NoTreatment" & treatment == 0 & fup < 5 
	replace wm_outcome = 0 if arm == "NoTreatment" & treatment == 0 & fup < 5

	** add 1 days to 0-survivors
	replace wm_fup= 1 if arm == "NoTreatment" & wm_fup==0

	order pseudo_key arm
	compress
	save "Sensitivity/SEN2/cloned/cloned `event' SEN2 bs`j'.dta", replace
}
}
*
* Split times
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "split" _col(70) "SEN2"
	}
	use "Sensitivity/SEN2/cloned/cloned `event' SEN2 bs`j'.dta", clear

	** times of event
	stset fup, fail(outcome) id(NewID)
	stsplit, at(failures)

	gen tstart = _t0
	gen tstop = _t

	** times of censoring
	gen TrialEmul_cens = 1-outcome
	stset fup, fail(TrialEmul_cens) id(NewID)
	stsplit, at(failures)

	replace tstart = _t0 if tstart<_t0 & _t0 != . & _t != .
	replace tstop = _t if tstop>_t & _t0 != . & _t != .

	order pseudo_key arm tstart tstop
	sort _all
	compress

	save "Sensitivity/SEN2/split/split `event' SEN2 bs`j'.dta", replace
}
}
*
* IPCW - Treatment arm
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "Treatment" _col(70) "SEN2"
	}
	use "Sensitivity/SEN2/split/split `event' SEN2 bs`j'.dta", clear
	keep if arm == "Treatment"
	replace outcome = . if TrialEmul_cens == .
	merge m:1 pseudo_key using "covid pregnant antiviral outpatient characteristics.dta", keep(3) nogen

	* adapt records to the long format
	sort NewID tstart
	bysort NewID (tstart): replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	capture stcox c.age#i.age_gp i.vaccine_status i.month_covid i.immuno_hist i.HELLP_bl i.GDM i.MMMI_bl i.preg_ht_preeclampsia_bl, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "Sensitivity/SEN2/Treatment/Treatment `event' SEN2 bs`j'.dta", replace
}
}
*
* IPCW - Control arm
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "NoTreatment" _col(70) "SEN2"
	}
	use "Sensitivity/SEN2/split/split `event' SEN2 bs`j'.dta", clear
	keep if arm == "NoTreatment"
	replace outcome = . if TrialEmul_cens == .
	merge m:1 pseudo_key using "covid pregnant antiviral outpatient characteristics.dta", keep(3) nogen

	* adapt records to the long format
	sort NewID tstart
	bysort NewID: replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	capture stcox c.age#i.age_gp i.vaccine_status i.month_covid i.immuno_hist i.HELLP_bl i.GDM i.MMMI_bl i.preg_ht_preeclampsia_bl, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "Sensitivity/SEN2/NoTreatment/Notreatment `event' SEN2 bs`j'.dta", replace
}
}
*
* Combine & Generate weights
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "Combine" _col(70) "SEN2"
	}
	use "Sensitivity/SEN2/Treatment/Treatment `event' SEN2 bs`j'.dta", clear
	append using "Sensitivity/SEN2/NoTreatment/Notreatment `event' SEN2 bs`j'.dta"

	// create a new ID variable for each clone in each arm
	tostring NewID, gen(NewID_str)
	replace NewID_str = "00000" + NewID_str if length(NewID_str)== 1
	replace NewID_str = "0000" + NewID_str if length(NewID_str)== 2
	replace NewID_str = "000" + NewID_str if length(NewID_str)== 3
	replace NewID_str = "00" + NewID_str if length(NewID_str)== 4
	replace NewID_str = "0" + NewID_str if length(NewID_str)== 5
	gen Anal_ID = "1" + NewID_str if arm == "Treatment" 
	replace Anal_ID = "2" + NewID_str if arm == "NoTreatment" 

	replace weight = 1 if wm_outcome == 1 & tstop == 5

	gen arm_value = 1 if arm == "Treatment"
	replace arm_value = 0 if arm_value == .

	keep pseudo_key tstart tstop fup_obs event time_to_treatment treatment bs outcome fup wm_fup wm_outcome TrialEmul_cens weight Anal_ID arm_value invalid
	destring Anal_ID, replace
	compress
	save "Sensitivity/SEN2/emulated/emulated `event' SEN2 bs`j'.dta", replace
}
}
*
* Generate weights
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "ipcw" _col(70) "SEN2"
	}
	use "Sensitivity/SEN2/emulated/emulated `event' SEN2 bs`j'.dta", clear
	rename weight _ipcw
	compress
	save "Sensitivity/SEN2/weight/weight `event' SEN2 bs`j'.dta", replace
}
}
*
cls
* Generate KM estimate
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "KM ipcw" _col(80) "SEN2"
	}
	use "Sensitivity/SEN2/weight/weight `event' SEN2 bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = _ipcw], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(fup bs invalid)
	save "Sensitivity/SEN2/KM ipcw/KM `event' SEN2 bs`j'.dta", replace
}
}
*
* Finalize bootstrap datasets
foreach event in covid_hosp MMMI c_section preterm_birth_comb {
clear
forvalues j = 0/500 {
	append using "Sensitivity/SEN2/KM ipcw/KM `event' SEN2 bs`j'.dta"
}
	gen hazard_s_w = 1 - KM_s_w
	gen hazard_ns_w = 1 - KM_ns_w
	gen odds_s_w = hazard_s_w/(1-hazard_s_w)
	gen odds_ns_w = hazard_ns_w/(1-hazard_ns_w)
	gen RR_w = hazard_s_w/hazard_ns_w
	gen diff_w = hazard_s_w - hazard_ns_w
	gen OR_w = odds_s_w / odds_ns_w
	compress
	save "Sensitivity/SEN2/KM ipcw `event' SEN2 bs_all.dta", replace
}
*
cls
* Cumulative incidence (N, %)
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	use "Sensitivity/SEN2/KM ipcw `event' SEN2 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum hazard_s_w if bs == 0, d
	scalar hazard_s_mean = r(mean)
	centile hazard_s_w if bs > 0, centile(2.5 97.5)
	scalar hazard_s_cil = r(c_1)
	scalar hazard_s_ciu = r(c_2)
	
	sum hazard_ns_w if bs == 0, d
	scalar hazard_ns_mean = r(mean)
	centile hazard_ns_w if bs > 0, centile(2.5 97.5)
	scalar hazard_ns_cil = r(c_1)
	scalar hazard_ns_ciu = r(c_2)
	
	noi di substr("`event'",1,13) _col(15) "SEN2" _col(30) hazard_s_mean _col(45) hazard_ns_mean
}
*
* Absolute risk reduction
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	use "Sensitivity/SEN2/KM ipcw `event' SEN2 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum diff_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile diff_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "SEN2" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
*
* Relative risk
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	use "Sensitivity/SEN2/KM ipcw `event' SEN2 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum RR_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile RR_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "SEN2" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
*

******************************************
*** Sen3 - Using cloning-censoring-IPCW method without prior matching by propensity-score

***Creates directory
capture mkdir "Sensitivity/SEN3"
capture mkdir "Sensitivity/SEN3/prepare"
capture mkdir "Sensitivity/SEN3/cloned"
capture mkdir "Sensitivity/SEN3/split"
capture mkdir "Sensitivity/SEN3/Treatment"
capture mkdir "Sensitivity/SEN3/NoTreatment"
capture mkdir "Sensitivity/SEN3/emulated"
capture mkdir "Sensitivity/SEN3/weight"
capture mkdir "Sensitivity/SEN3/KM ipcw"
*

*** Prepare dataset for analysis
cls
qui foreach event in covid_hosp MMMI {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {	
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "setup" _col(70) "SEN3"
	}
	* Setup for trial emulation
	use "covid pregnant antiviral outpatient.dta", clear
	keep if group < .

	* events
	merge 1:1 pseudo_key using "date_covid_hosp", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_MMMI", keep(1 3) nogen
	
	if `j' > 0 {
		set seed `j'
		bsample
	}
	replace date_paxlovid = . if date_paxlovid >= date_admission
	* organize / rename / generate variables
	gen date_last_fu = min(date_death, mdy(02,12,2023), date_baseline + 28)
	gen date_event = date_`event'
	gen event = inrange(date_event, date_baseline, date_last_fu)
	gen fup_obs = min(date_event-date_baseline, 28) if event == 1
	replace fup_obs = min(date_last_fu-date_baseline, 28) if event == 0
	tab fup_obs
	gen time_to_treatment = date_paxlovid - date_baseline
	tab time_to_treatment
	gen treatment = inrange(date_paxlovid - date_baseline, 0, 5)
	* keep necessary variables
	keep pseudo_key fup_obs event time_to_treatment treatment
	gen bs = `j'
	compress
	save "Sensitivity/SEN3/prepare/prepare `event' SEN3 bs`j'.dta", replace
}
}
*
qui foreach event in c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "setup" _col(70) "SEN3"
	}
	* Setup for trial emulation
	use "covid pregnant antiviral outpatient.dta", clear
	keep if group < .

	* events
	merge 1:1 pseudo_key using "date_c_section", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_preterm_birth", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_stillbirth", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_neonatal_death", keep(1 3) nogen
	
	if `j' > 0 {
		set seed `j'
		bsample
	}
	replace date_paxlovid = . if date_paxlovid >= date_admission
	* organize / rename / generate variables
	gen date_last_fu = min(date_death, mdy(02,12,2023))
	gen date_event = date_`event'
	gen event = date_event < .
	gen fup_obs = date_event-date_baseline if event == 1
	replace fup_obs = date_last_fu-date_baseline if event == 0
	tab fup_obs
	gen time_to_treatment = date_paxlovid - date_baseline
	tab time_to_treatment
	gen treatment = inrange(date_paxlovid - date_baseline, 0, 5)
	* keep necessary variables
	keep pseudo_key fup_obs event time_to_treatment treatment
	gen bs = `j'
	compress
	save "Sensitivity/SEN3/prepare/prepare `event' SEN3 bs`j'.dta", replace
}
}
*
cls
***** Bootstrap
*** Cloning & censoring
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "cloning" _col(70) "SEN3"
	}
	* Prepare dataset for analysis
	use "Sensitivity/SEN3/prepare/prepare `event' SEN3 bs`j'.dta", replace
	stset fup_obs, failure(event)

	* Arm A: no treatment within 5 days (control: non-exposed group)
	gen outcomeA = _d // _d = `event'
	gen fupA = _t // _t = follow up time

	/// if the patient received treatment within 5 days:
	/// 1. no event outcome, since the patient survived till censoring (treatment)
	replace outcomeA = 0 if treatment==1 & time_to_treatment <=5 
	/// 2. follow up is censored at treatment
	replace fupA = time_to_treatment if treatment==1 & time_to_treatment <=5

	* Arm B: treatment within 5 days (treated: exposed group)
	gen outcomeB = _d 
	gen fupB = _t 

	/// if the patient survived the first 5 days and did not receive treatment within 5 days:
	/// 1. no event outcome if the patient survived the first 5 days
	replace outcomeB = 0 if (treatment==0 & _t>5) | (treatment==1 & time_to_treatment >5 & time_to_treatment !=.)
	/// 2. follow up is censored at 5 days
	replace fupB = 5 if (treatment==0 & _t>5) | (treatment==1 & time_to_treatment >5 & time_to_treatment != .)

	** append clones 
	preserve
		drop outcomeB fupB
		rename outcomeA outcome
		rename fupA fup
		gen arm = "NoTreatment"
		tempfile a
		save "`a'", replace
	restore
		drop outcomeA fupA
		rename outcomeB outcome
		rename fupB fup
		gen arm = "Treatment"	
		cap append using "`a'"

	// Weight models

	sort _all
	gen NewID = _n

	** add 1 day to 0-survivors
	replace fup= 1 if fup==0

	** Weight model: define survival time and event indicator	
	* treatment Arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	gen wm_fup = time_to_treatment if arm == "Treatment" & time_to_treatment<=5 & time_to_treatment!=. & treatment == 1 
	gen wm_outcome = 0 if arm == "Treatment" & time_to_treatment<=5 & time_to_treatment!=. & treatment == 1 

	** Case 2: they deviate at 5 days
	replace wm_fup = 5 if arm == "Treatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1))
	replace wm_outcome = 1 if arm == "Treatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "Treatment" & treatment == 0 & fup < 5
	replace wm_outcome = 0 if arm == "Treatment" & treatment == 0 & fup < 5
	** add 1 days to 0-survivors
	replace wm_fup= 1 if arm == "Treatment" & wm_fup==0 

	* No treatment Arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	replace wm_fup = time_to_treatment if arm == "NoTreatment" & time_to_treatment<=5 & treatment == 1 
	replace wm_outcome = 1 if arm == "NoTreatment" & time_to_treatment<=5 & treatment == 1 

	** Case 2: they deviate at 5 days
	replace wm_fup = 5 if arm == "NoTreatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 
	replace wm_outcome = 0 if arm == "NoTreatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "NoTreatment" & treatment == 0 & fup < 5 
	replace wm_outcome = 0 if arm == "NoTreatment" & treatment == 0 & fup < 5

	** add 1 days to 0-survivors
	replace wm_fup= 1 if arm == "NoTreatment" & wm_fup==0

	order pseudo_key arm
	compress
	save "Sensitivity/SEN3/cloned/cloned `event' SEN3 bs`j'.dta", replace
}
}
*
* Split times
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "split" _col(70) "SEN3"
	}
	use "Sensitivity/SEN3/cloned/cloned `event' SEN3 bs`j'.dta", clear

	** times of event
	stset fup, fail(outcome) id(NewID)
	stsplit, at(failures)

	gen tstart = _t0
	gen tstop = _t

	** times of censoring
	gen TrialEmul_cens = 1-outcome
	stset fup, fail(TrialEmul_cens) id(NewID)
	stsplit, at(failures)

	replace tstart = _t0 if tstart<_t0 & _t0 != . & _t != .
	replace tstop = _t if tstop>_t & _t0 != . & _t != .

	order pseudo_key arm tstart tstop
	sort _all
	compress

	save "Sensitivity/SEN3/split/split `event' SEN3 bs`j'.dta", replace
}
}
*
* IPCW - Treatment arm
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "Treatment" _col(70) "SEN3"
	}
	use "Sensitivity/SEN3/split/split `event' SEN3 bs`j'.dta", clear
	keep if arm == "Treatment"
	replace outcome = . if TrialEmul_cens == .
	merge m:1 pseudo_key using "covid pregnant antiviral outpatient characteristics.dta", keep(3) nogen

	* adapt records to the long format
	sort NewID tstart
	bysort NewID (tstart): replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	capture stcox c.age#i.age_gp i.vaccine_status i.month_covid i.immuno_hist i.HELLP_bl i.GDM i.MMMI_bl i.preg_ht_preeclampsia_bl, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "Sensitivity/SEN3/Treatment/Treatment `event' SEN3 bs`j'.dta", replace
}
}
*
* IPCW - Control arm
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "NoTreatment" _col(70) "SEN3"
	}
	use "Sensitivity/SEN3/split/split `event' SEN3 bs`j'.dta", clear
	keep if arm == "NoTreatment"
	replace outcome = . if TrialEmul_cens == .
	merge m:1 pseudo_key using "covid pregnant antiviral outpatient characteristics.dta", keep(3) nogen

	* adapt records to the long format
	sort NewID tstart
	bysort NewID: replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	capture stcox c.age#i.age_gp i.vaccine_status i.month_covid i.immuno_hist i.HELLP_bl i.GDM i.MMMI_bl i.preg_ht_preeclampsia_bl, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "Sensitivity/SEN3/NoTreatment/Notreatment `event' SEN3 bs`j'.dta", replace
}
}
*
* Combine & Generate weights
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "Combine" _col(70) "SEN3"
	}
	use "Sensitivity/SEN3/Treatment/Treatment `event' SEN3 bs`j'.dta", clear
	append using "Sensitivity/SEN3/NoTreatment/Notreatment `event' SEN3 bs`j'.dta"

	// create a new ID variable for each clone in each arm
	tostring NewID, gen(NewID_str)
	replace NewID_str = "00000" + NewID_str if length(NewID_str)== 1
	replace NewID_str = "0000" + NewID_str if length(NewID_str)== 2
	replace NewID_str = "000" + NewID_str if length(NewID_str)== 3
	replace NewID_str = "00" + NewID_str if length(NewID_str)== 4
	replace NewID_str = "0" + NewID_str if length(NewID_str)== 5
	gen Anal_ID = "1" + NewID_str if arm == "Treatment" 
	replace Anal_ID = "2" + NewID_str if arm == "NoTreatment" 

	replace weight = 1 if wm_outcome == 1 & tstop == 5

	gen arm_value = 1 if arm == "Treatment"
	replace arm_value = 0 if arm_value == .

	keep pseudo_key tstart tstop fup_obs event time_to_treatment treatment bs outcome fup wm_fup wm_outcome TrialEmul_cens weight Anal_ID arm_value invalid
	destring Anal_ID, replace

	compress
	save "Sensitivity/SEN3/emulated/emulated `event' SEN3 bs`j'.dta", replace
}
}
*
* Generate weights
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "ipcw" _col(70) "SEN3"
	}
	use "Sensitivity/SEN3/emulated/emulated `event' SEN3 bs`j'.dta", clear
	rename weight _ipcw
	compress
	save "Sensitivity/SEN3/weight/weight `event' SEN3 bs`j'.dta", replace
}
}
*
cls
* Generate KM estimate
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "KM ipcw" _col(80) "SEN3"
	}
	use "Sensitivity/SEN3/weight/weight `event' SEN3 bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = _ipcw], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(fup bs invalid)
	save "Sensitivity/SEN3/KM ipcw/KM `event' SEN3 bs`j'.dta", replace
}
}
*
* Finalize bootstrap datasets
foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	clear
forvalues j = 0/500 {
	append using "Sensitivity/SEN3/KM ipcw/KM `event' SEN3 bs`j'.dta"
}
	gen hazard_s_w = 1 - KM_s_w
	gen hazard_ns_w = 1 - KM_ns_w
	gen odds_s_w = hazard_s_w/(1-hazard_s_w)
	gen odds_ns_w = hazard_ns_w/(1-hazard_ns_w)
	gen RR_w = hazard_s_w/hazard_ns_w
	gen diff_w = hazard_s_w - hazard_ns_w
	gen OR_w = odds_s_w / odds_ns_w
	compress
	save "Sensitivity/SEN3/KM ipcw `event' SEN3 bs_all.dta", replace
}
*
cls
* Cumulative incidence (N, %)
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	use "Sensitivity/SEN3/KM ipcw `event' SEN3 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum hazard_s_w if bs == 0, d
	scalar hazard_s_mean = r(mean)
	centile hazard_s_w if bs > 0, centile(2.5 97.5)
	scalar hazard_s_cil = r(c_1)
	scalar hazard_s_ciu = r(c_2)
	
	sum hazard_ns_w if bs == 0, d
	scalar hazard_ns_mean = r(mean)
	centile hazard_ns_w if bs > 0, centile(2.5 97.5)
	scalar hazard_ns_cil = r(c_1)
	scalar hazard_ns_ciu = r(c_2)
	
	noi di substr("`event'",1,13) _col(15) "SEN3" _col(30) hazard_s_mean _col(45) hazard_ns_mean
}
*
* Absolute risk reduction
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	use "Sensitivity/SEN3/KM ipcw `event' SEN3 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum diff_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile diff_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "SEN3" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
*
* Relative risk
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	use "Sensitivity/SEN3/KM ipcw `event' SEN3 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum RR_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile RR_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "SEN3" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
*

******************************************
*** Sen4 - Performing propensity-score matching by day of nirmatrelvir/ritonavir initiation

***Creates directory
capture mkdir "Sensitivity/SEN4"
capture mkdir "Sensitivity/SEN4/prepare"
capture mkdir "Sensitivity/SEN4/cloned"
capture mkdir "Sensitivity/SEN4/split"
capture mkdir "Sensitivity/SEN4/Treatment"
capture mkdir "Sensitivity/SEN4/NoTreatment"
capture mkdir "Sensitivity/SEN4/emulated"
capture mkdir "Sensitivity/SEN4/weight"
capture mkdir "Sensitivity/SEN4/KM ipcw"
*

*** Define SEN4
forvalues i = 1/1 {
	use "covid_pregnant_patient_list.dta", clear
	merge 1:1 pseudo_key using "eGFR_30_hist.dta", keep(1 3) nogen
	merge 1:1 pseudo_key using "renal_liver_hist.dta", keep(1 3) nogen
	merge 1:1 pseudo_key using "case.dta", keep(1 3) keepusing(age_case asymptomatic date_onset recovered imported) nogen

	drop if asymptomatic == 1

	gen elderly_home = elderly_home_ind == "Y"

	drop if (paxlovid_DDI_3 == 1 | paxlovid_DDI_2 == 1)
	drop if (eGFR_30_hist == 1 | dialysis_hist == 1 | renal_trans_hist == 1)
	drop if cirrhosis_hist == 1 | liver_cancer_hist == 1 | liver_trans_hist == 1
	drop if age <= 40 & elderly_home == 1
	keep if elderly_home == 0

	rename gp_P group
	tab group

	gen day_trt_covid = date_paxlovid - date_baseline if group == 1
	tab day_trt_covid group

	gen covid_infect = (date_covid >= mdy(9,1,2022))

	save "Sensitivity/SEN4/covid pregnant antiviral outpatient SEN4.dta", replace

	* Day 0-5
	qui forvalues k = 0/5 {
		use "Sensitivity/SEN4/covid pregnant antiviral outpatient SEN4.dta", clear
		if `k' == 1 {
			merge 1:1 pseudo_key using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day0.dta", keep(1) keepusing(pseudo_key day_trt_covid) nogen
		}
		if `k' == 2 {
			merge 1:1 pseudo_key using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day0.dta", keep(1) keepusing(pseudo_key day_trt_covid) nogen
			merge 1:1 pseudo_key using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day1.dta", keep(1) keepusing(pseudo_key day_trt_covid) nogen
		}
		if `k' == 3 {
			merge 1:1 pseudo_key using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day0.dta", keep(1) keepusing(pseudo_key day_trt_covid) nogen
			merge 1:1 pseudo_key using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day1.dta", keep(1) keepusing(pseudo_key day_trt_covid) nogen
			merge 1:1 pseudo_key using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day2.dta", keep(1) keepusing(pseudo_key day_trt_covid) nogen
		}
		if `k' == 4 {
			merge 1:1 pseudo_key using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day0.dta", keep(1) keepusing(pseudo_key day_trt_covid) nogen
			merge 1:1 pseudo_key using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day1.dta", keep(1) keepusing(pseudo_key day_trt_covid) nogen
			merge 1:1 pseudo_key using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day2.dta", keep(1) keepusing(pseudo_key day_trt_covid) nogen
			merge 1:1 pseudo_key using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day3.dta", keep(1) keepusing(pseudo_key day_trt_covid) nogen
		}
		if `k' == 5 {
			merge 1:1 pseudo_key using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day0.dta", keep(1) keepusing(pseudo_key day_trt_covid) nogen
			merge 1:1 pseudo_key using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day1.dta", keep(1) keepusing(pseudo_key day_trt_covid) nogen
			merge 1:1 pseudo_key using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day2.dta", keep(1) keepusing(pseudo_key day_trt_covid) nogen
			merge 1:1 pseudo_key using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day3.dta", keep(1) keepusing(pseudo_key day_trt_covid) nogen
			merge 1:1 pseudo_key using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day4.dta", keep(1) keepusing(pseudo_key day_trt_covid) nogen
		}
		merge 1:1 pseudo_key using "covid pregnant antiviral outpatient characteristics.dta", keep(1 3) nogen
		keep if day_trt_covid == `k' | group == 0

		local mi = "logit group age#i.age_gp c.date_covid#i.covid_infect i.vaccine_status i.immuno_hist i.HELLP_bl i.GDM i.MMMI_bl i.preg_ht_preeclampsia_bl"

		`mi'
		predict prob_group
		gen xb_group = ln(prob_group / (1-prob_group))

		rename group group_old

		noi di "matching 1:10"
		set seed 123456
		calipmatch, gen(group_10_match) case(group_old) max(10) calipermatch(xb_group) caliperwidth(.05)
		bysort pseudo_key (group_10_match): replace group_10_match = group_10_match[1]
		gen group_10 = group_old if group_10_match < .

		keep if group_10 < .

		noi di "day_trt = `k'"
		noi tab group_10
		
		compress
		save "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day`k'.dta", replace
	}
	*

	use "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day0.dta", clear
	append using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day1.dta"
	append using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day2.dta"
	append using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day3.dta"
	append using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day4.dta"
	append using "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4_day5.dta"

	tab group_10

	save "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4 dayk_combined.dta", replace
}
*
*** Prepare dataset for analysis
qui foreach event in covid_hosp MMMI {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "setup" _col(70) "SEN4"
	}
	* Setup for trial emulation
	use "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4 dayk_combined.dta", clear
	rename group_10 group
	keep if group < .
	
	* events
	merge 1:1 pseudo_key using "date_covid_hosp", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_MMMI", keep(1 3) nogen
	
	if `j' > 0 {
		set seed `j'
		bsample
	}
	replace date_paxlovid = . if date_paxlovid >= date_admission
	* organize / rename / generate variables
	gen date_last_fu = min(date_death, mdy(02,12,2023), date_baseline + 28)
	gen date_event = date_`event'
	gen event = inrange(date_event, date_baseline, date_last_fu)
	gen fup_obs = min(date_event-date_baseline, 28) if event == 1
	replace fup_obs = min(date_last_fu-date_baseline, 28) if event == 0
	tab fup_obs
	gen time_to_treatment = date_paxlovid - date_baseline
	tab time_to_treatment
	gen treatment = inrange(date_paxlovid - date_baseline, 0, 5)
	* keep necessary variables
	keep pseudo_key fup_obs event time_to_treatment treatment
	gen bs = `j'
	compress
	save "Sensitivity/SEN4/prepare/prepare `event' SEN4 bs`j'.dta", replace
}
}
*
qui foreach event in c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "setup" _col(70) "SEN4"
	}
	* Setup for trial emulation
	use "Sensitivity/SEN4/covid pregnant antiviral outpatient matched SEN4 dayk_combined.dta", clear
	rename group_10 group
	keep if group < .
	
	* events
	merge 1:1 pseudo_key using "date_c_section", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_preterm_birth", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_stillbirth", keep(1 3) nogen
	merge 1:1 pseudo_key using "date_neonatal_death", keep(1 3) nogen
	
	if `j' > 0 {
		set seed `j'
		bsample
	}
	replace date_paxlovid = . if date_paxlovid >= date_admission
	* organize / rename / generate variables
	gen date_last_fu = min(date_death, mdy(02,12,2023))
	gen date_event = date_`event'
	gen event = inrange(date_event, date_baseline, date_last_fu)
	gen fup_obs = date_event-date_baseline if event == 1
	replace fup_obs = date_last_fu-date_baseline if event == 0
	tab fup_obs
	gen time_to_treatment = date_paxlovid - date_baseline
	tab time_to_treatment
	gen treatment = inrange(date_paxlovid - date_baseline, 0, 5)
	* keep necessary variables
	keep pseudo_key fup_obs event time_to_treatment treatment
	gen bs = `j'
	compress
	save "Sensitivity/SEN4/prepare/prepare `event' SEN4 bs`j'.dta", replace
}
}
*
cls
***** Bootstrap
*** Cloning & censoring
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "cloning" _col(70) "SEN4" 
	}
	* Prepare dataset for analysis
	use "Sensitivity/SEN4/prepare/prepare `event' SEN4 bs`j'.dta", replace
	stset fup_obs, failure(event)

	* Arm A: no treatment within 5 days (control: non-exposed group)
	gen outcomeA = _d // _d = `event'
	gen fupA = _t // _t = follow up time

	/// if the patient received treatment within 5 days:
	/// 1. no event outcome, since the patient survived till censoring (treatment)
	replace outcomeA = 0 if treatment==1 & time_to_treatment <=5 
	/// 2. follow up is censored at treatment
	replace fupA = time_to_treatment if treatment==1 & time_to_treatment <=5

	* Arm B: treatment within 5 days (treated: exposed group)
	gen outcomeB = _d 
	gen fupB = _t 

	/// if the patient survived the first 5 days and did not receive treatment within 5 days:
	/// 1. no event outcome if the patient survived the first 5 days
	replace outcomeB = 0 if (treatment==0 & _t>5) | (treatment==1 & time_to_treatment >5 & time_to_treatment !=.)
	/// 2. follow up is censored at 5 days
	replace fupB = 5 if (treatment==0 & _t>5) | (treatment==1 & time_to_treatment >5 & time_to_treatment != .)

	** append clones 
	preserve
		drop outcomeB fupB
		rename outcomeA outcome
		rename fupA fup
		gen arm = "NoTreatment"
		tempfile a
		save "`a'", replace
	restore
		drop outcomeA fupA
		rename outcomeB outcome
		rename fupB fup
		gen arm = "Treatment"	
		cap append using "`a'"

	// Weight models

	sort _all
	gen NewID = _n

	** add 1 day to 0-survivors
	replace fup= 1 if fup==0

	** Weight model: define survival time and event indicator	
	* treatment Arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	gen wm_fup = time_to_treatment if arm == "Treatment" & time_to_treatment<=5 & time_to_treatment!=. & treatment == 1 
	gen wm_outcome = 0 if arm == "Treatment" & time_to_treatment<=5 & time_to_treatment!=. & treatment == 1 

	** Case 2: they deviate at 5 days
	replace wm_fup = 5 if arm == "Treatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1))
	replace wm_outcome = 1 if arm == "Treatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "Treatment" & treatment == 0 & fup < 5
	replace wm_outcome = 0 if arm == "Treatment" & treatment == 0 & fup < 5
	** add 1 days to 0-survivors
	replace wm_fup= 1 if arm == "Treatment" & wm_fup==0 

	* No treatment Arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	replace wm_fup = time_to_treatment if arm == "NoTreatment" & time_to_treatment<=5 & treatment == 1 
	replace wm_outcome = 1 if arm == "NoTreatment" & time_to_treatment<=5 & treatment == 1 

	** Case 2: they deviate at 5 days
	replace wm_fup = 5 if arm == "NoTreatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 
	replace wm_outcome = 0 if arm == "NoTreatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "NoTreatment" & treatment == 0 & fup < 5 
	replace wm_outcome = 0 if arm == "NoTreatment" & treatment == 0 & fup < 5

	** add 1 days to 0-survivors
	replace wm_fup= 1 if arm == "NoTreatment" & wm_fup==0

	order pseudo_key arm
	compress
	save "Sensitivity/SEN4/cloned/cloned `event' SEN4 bs`j'.dta", replace
}
}
*
* Split times
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "split" _col(70) "SEN4"
	}
	use "Sensitivity/SEN4/cloned/cloned `event' SEN4 bs`j'.dta", clear

	** times of event
	stset fup, fail(outcome) id(NewID)
	stsplit, at(failures)

	gen tstart = _t0
	gen tstop = _t

	** times of censoring
	gen TrialEmul_cens = 1-outcome
	stset fup, fail(TrialEmul_cens) id(NewID)
	stsplit, at(failures)

	replace tstart = _t0 if tstart<_t0 & _t0 != . & _t != .
	replace tstop = _t if tstop>_t & _t0 != . & _t != .

	order pseudo_key arm tstart tstop
	sort _all
	compress

	save "Sensitivity/SEN4/split/split `event' SEN4 bs`j'.dta", replace
}
}
*
* IPCW - Treatment arm
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "Treatment" _col(70) "SEN4"
	}
	use "Sensitivity/SEN4/split/split `event' SEN4 bs`j'.dta", clear
	keep if arm == "Treatment"
	replace outcome = . if TrialEmul_cens == .
	merge m:1 pseudo_key using "covid pregnant antiviral outpatient characteristics.dta", keep(3) nogen

	* adapt records to the long format
	sort NewID tstart
	bysort NewID (tstart): replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	capture stcox c.age#i.age_gp i.vaccine_status i.month_covid i.immuno_hist i.HELLP_bl i.GDM i.MMMI_bl i.preg_ht_preeclampsia_bl, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "Sensitivity/SEN4/Treatment/Treatment `event' SEN4 bs`j'.dta", replace
}
}
*
* IPCW - Control arm
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "NoTreatment" _col(70) "SEN4"
	}
	use "Sensitivity/SEN4/split/split `event' SEN4 bs`j'.dta", clear
	keep if arm == "NoTreatment"
	replace outcome = . if TrialEmul_cens == .
	merge m:1 pseudo_key using "covid pregnant antiviral outpatient characteristics.dta", keep(3) nogen

	* adapt records to the long format
	sort NewID tstart
	bysort NewID: replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	capture stcox c.age#i.age_gp i.vaccine_status i.month_covid i.immuno_hist i.HELLP_bl i.GDM i.MMMI_bl i.preg_ht_preeclampsia_bl, efron
	if _rc == 0 {
		predict ch2_2, basech
		gen s0_2=exp(-ch2_2)
		predict xbCox2_2, xb
		gen weight = 1/s0_2^(exp(xbCox2_2)) 
		drop _st-xbCox2_2 
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen weight = .
		noi di "invalid"
	}
	keep pseudo_key arm tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "Sensitivity/SEN4/NoTreatment/Notreatment `event' SEN4 bs`j'.dta", replace
}
}
*
* Combine & Generate weights
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "Combine" _col(70) "SEN4"
	}
	use "Sensitivity/SEN4/Treatment/Treatment `event' SEN4 bs`j'.dta", clear
	append using "Sensitivity/SEN4/NoTreatment/Notreatment `event' SEN4 bs`j'.dta"

	// create a new ID variable for each clone in each arm
	tostring NewID, gen(NewID_str)
	replace NewID_str = "00000" + NewID_str if length(NewID_str)== 1
	replace NewID_str = "0000" + NewID_str if length(NewID_str)== 2
	replace NewID_str = "000" + NewID_str if length(NewID_str)== 3
	replace NewID_str = "00" + NewID_str if length(NewID_str)== 4
	replace NewID_str = "0" + NewID_str if length(NewID_str)== 5
	gen Anal_ID = "1" + NewID_str if arm == "Treatment" 
	replace Anal_ID = "2" + NewID_str if arm == "NoTreatment" 

	replace weight = 1 if wm_outcome == 1 & tstop == 5

	gen arm_value = 1 if arm == "Treatment"
	replace arm_value = 0 if arm_value == .

	keep pseudo_key tstart tstop fup_obs event time_to_treatment treatment bs outcome fup wm_fup wm_outcome TrialEmul_cens weight Anal_ID arm_value invalid
	destring Anal_ID, replace
	compress
	save "Sensitivity/SEN4/emulated/emulated `event' SEN4 bs`j'.dta", replace
}
}
*
* Generate weights
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "ipcw" _col(70) "SEN4"
	}
	use "Sensitivity/SEN4/emulated/emulated `event' SEN4 bs`j'.dta", clear
	rename weight _ipcw
	compress
	save "Sensitivity/SEN4/weight/weight `event' SEN4 bs`j'.dta", replace
}
}
*
cls
* Generate KM estimate
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "KM ipcw" _col(80) "SEN4"
	}
	use "Sensitivity/SEN4/weight/weight `event' SEN4 bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = _ipcw], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(fup bs invalid)
	save "Sensitivity/SEN4/KM ipcw/KM `event' SEN4 bs`j'.dta", replace
}
}
*
* Finalize bootstrap datasets
foreach event in covid_hosp MMMI c_section preterm_birth_comb {
clear
forvalues j = 0/500 {
	append using "Sensitivity/SEN4/KM ipcw/KM `event' SEN4 bs`j'.dta"
}
	gen hazard_s_w = 1 - KM_s_w
	gen hazard_ns_w = 1 - KM_ns_w
	gen odds_s_w = hazard_s_w/(1-hazard_s_w)
	gen odds_ns_w = hazard_ns_w/(1-hazard_ns_w)
	gen RR_w = hazard_s_w/hazard_ns_w
	gen diff_w = hazard_s_w - hazard_ns_w
	gen OR_w = odds_s_w / odds_ns_w
	compress
	save "Sensitivity/SEN4/KM ipcw `event' SEN4 bs_all.dta", replace
}
*
cls
* Cumulative incidence (N, %)
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	use "Sensitivity/SEN4/KM ipcw `event' SEN4 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum hazard_s_w if bs == 0, d
	scalar hazard_s_mean = r(mean)
	centile hazard_s_w if bs > 0, centile(2.5 97.5)
	scalar hazard_s_cil = r(c_1)
	scalar hazard_s_ciu = r(c_2)
	
	sum hazard_ns_w if bs == 0, d
	scalar hazard_ns_mean = r(mean)
	centile hazard_ns_w if bs > 0, centile(2.5 97.5)
	scalar hazard_ns_cil = r(c_1)
	scalar hazard_ns_ciu = r(c_2)
	
	noi di substr("`event'",1,13) _col(15) "SEN4" _col(30) hazard_s_mean _col(45) hazard_ns_mean
}
*
* Absolute risk reduction
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	use "Sensitivity/SEN4/KM ipcw `event' SEN4 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum diff_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile diff_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "SEN4" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
*
* Relative risk
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	use "Sensitivity/SEN4/KM ipcw `event' SEN4 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum RR_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile RR_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "SEN4" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
*

******************************************
*** Sen5 - Using pooled logistic regression for IPCW model
capture mkdir "Sensitivity/SEN5"
capture mkdir "Sensitivity/SEN5/Treatment"
capture mkdir "Sensitivity/SEN5/Treatment pr_censor"
capture mkdir "Sensitivity/SEN5/NoTreatment"
capture mkdir "Sensitivity/SEN5/NoTreatment pr_censor"
capture mkdir "Sensitivity/SEN5/emulated"
capture mkdir "Sensitivity/SEN5/weight"
capture mkdir "Sensitivity/SEN5/KM ipcw"
*

* IPCW - Treatment arm
* Estimate proablility of censor
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	noi di "`event'" _col(15) "subgp1" _col(30) "bs`j'" _col(45) "Treatment pr_censor" _col(70) "SEN5"
	use "matched_10/subgp1/split/split `event' subgp_1 bs`j'.dta", clear
	keep if arm == "Treatment"

	* adapt records to the long format
	sort NewID tstart
	bysort NewID (tstart): replace wm_outcome = 0 if _n!=_N
	keep if tstart == 4

	merge m:1 pseudo_key using "covid pregnant antiviral outpatient characteristics.dta", keep(3) nogen

	* Weight model:
	capture logit wm_outcome c.age#i.age_gp i.vaccine_status i.month_covid i.immuno_hist i.HELLP_bl i.GDM i.MMMI_bl i.preg_ht_preeclampsia_bl
		if _rc == 0 {
			predict pr_censor if e(sample)
			gen invalid = 0
		}
		else {
			gen invalid = 1
			gen pr_censor = .
			noi di "invalid"
		}
	keep pseudo_key arm tstart pr_censor invalid bs
	bysort _all : keep if _n == 1
	compress
	save "Sensitivity/SEN5/Treatment pr_censor/Treatment pr_censor `event' SEN5 bs`j'.dta", replace
}
}
*
* Combine with daily dataset
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	noi di "subgp_1" _col(15) "bs`j'" _col(30) "Treatment" _col(55) "SEN5"
	use "matched_10/subgp1/split/split `event' subgp_1 bs`j'.dta", clear
	keep if arm == "Treatment"
	
	* fill in prob censor
	merge m:1 pseudo_key arm tstart using "Sensitivity/SEN5/Treatment pr_censor/Treatment pr_censor `event' SEN5 bs`j'.dta", keep(1 3) keepusing(pr_censor invalid) nogen
	replace pr_censor = 0 if pr_censor == .
	egen invalid_max = max(invalid)
	replace invalid = invalid_max
	drop invalid_max

	* expand to full time
	bysort pseudo_key tstart : gen n = _n
	gen diff = tstop - tstart
	expand diff
	sort pseudo_key n tstart
	bysort pseudo_key n tstart : replace tstart = tstart + _n - 1
	replace tstop = tstart + 1
	drop diff n
	
	* generate probability of uncensored
	gen pr_uncensor = 1 - pr_censor
	* generate denominator - cummulative probability of remained uncensored
	gen censdenom = pr_uncensor
	sort pseudo_key tstart
	by pseudo_key: replace censdenom = censdenom * censdenom[_n-1] if _n!=1

	keep pseudo_key arm tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens pr_censor pr_uncensor censdenom invalid
	compress
	save "Sensitivity/SEN5/Treatment/Treatment `event' SEN5 bs`j'.dta", replace
}
}
*
* IPCW - Control arm
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "NoTreatment pr_censor" _col(80) "SEN5"
	use "matched_10/subgp1/split/split `event' subgp_1 bs`j'.dta", clear
	keep if arm == "NoTreatment"
	replace outcome = . if TrialEmul_cens == .
	merge m:1 pseudo_key using "covid pregnant antiviral outpatient characteristics.dta", keep(3) nogen

	* expand to full time
	bysort pseudo_key tstart : gen n = _n
	gen diff = tstop - tstart
	expand diff
	sort pseudo_key n tstart
	bysort pseudo_key n tstart : replace tstart = tstart + _n - 1
	replace tstop = tstart + 1
	drop diff n

	* adapt records to the long format
	sort NewID tstart
	bysort NewID: replace wm_outcome = 0 if _n!=_N

	* Weight model:
	capture xi: logistic wm_outcome c.age#i.age_gp i.vaccine_status i.month_covid i.immuno_hist i.HELLP_bl i.GDM i.MMMI_bl i.preg_ht_preeclampsia_bl
	if _rc == 0 {
		predict pr_censor if e(sample)
		gen invalid = 0
	}
	else {
		gen invalid = 1
		gen pr_censor = .
		noi di "invalid"
	}
	keep pseudo_key arm tstart pr_censor invalid bs
	bysort _all : keep if _n == 1
	compress
	save "Sensitivity/SEN5/NoTreatment pr_censor/Notreatment pr_censor `event' SEN5 bs`j'.dta", replace
	}
}
*
* Combine with daily dataset
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	noi di "subgp_1" _col(15) "bs`j'" _col(30) "NoTreatment" _col(55) "SEN5"
	use "matched_10/subgp1/split/split `event' subgp_1 bs`j'.dta", clear
	keep if arm == "NoTreatment"
	
	* fill in prob censor
	merge m:1 pseudo_key arm tstart using "Sensitivity/SEN5/NoTreatment pr_censor/NoTreatment pr_censor `event' SEN5 bs`j'.dta", keep(1 3) keepusing(pr_censor invalid) nogen
	replace pr_censor = 0 if pr_censor == .
	egen invalid_max = max(invalid)
	replace invalid = invalid_max
	drop invalid_max

	* expand to full time
	bysort pseudo_key tstart : gen n = _n
	gen diff = tstop - tstart
	expand diff
	sort pseudo_key n tstart
	bysort pseudo_key n tstart : replace tstart = tstart + _n - 1
	replace tstop = tstart + 1
	drop diff n
	
	* generate probability of uncensored
	gen pr_uncensor = 1 - pr_censor
	* generate denominator - cummulative probability of remained uncensored
	gen censdenom = pr_uncensor
	sort pseudo_key tstart
	by pseudo_key: replace censdenom = censdenom * censdenom[_n-1] if _n!=1

	keep pseudo_key arm tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens pr_censor pr_uncensor censdenom invalid
	compress
	save "Sensitivity/SEN5/NoTreatment/NoTreatment `event' SEN5 bs`j'.dta", replace
}
}
*
* Combine & Generate weights
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	noi di "subgp_1" _col(15) "bs`j'" _col(30) "Combine" _col(55) "SEN5"
	use "Sensitivity/SEN5/Treatment/Treatment `event' SEN5 bs`j'.dta", clear
	append using "Sensitivity/SEN5/NoTreatment/Notreatment `event' SEN5 bs`j'.dta"

	// create a new ID variable for each clone in each arm
	tostring NewID, gen(NewID_str)
	replace NewID_str = "00000" + NewID_str if length(NewID_str)== 1
	replace NewID_str = "0000" + NewID_str if length(NewID_str)== 2
	replace NewID_str = "000" + NewID_str if length(NewID_str)== 3
	replace NewID_str = "00" + NewID_str if length(NewID_str)== 4
	replace NewID_str = "0" + NewID_str if length(NewID_str)== 5
	gen Anal_ID = "1" + NewID_str if arm == "Treatment" 
	replace Anal_ID = "2" + NewID_str if arm == "NoTreatment" 

	gen arm_value = 1 if arm == "Treatment"
	replace arm_value = 0 if arm_value == .

	keep pseudo_key tstart tstop fup_obs event time_to_treatment treatment bs outcome fup wm_fup wm_outcome TrialEmul_cens pr_censor pr_uncensor censdenom Anal_ID arm_value invalid
	destring Anal_ID, replace

	compress
	save "Sensitivity/SEN5/emulated/emulated `event' SEN5 bs`j'.dta", replace
}
}
*
* Generate weights
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "ipcw" _col(80) "SEN5"
	}
	use "Sensitivity/SEN5/emulated/emulated `event' SEN5 bs`j'.dta", clear
	bysort arm_value tstart : gen N_group = _N
	bysort arm_value tstart : egen N_censored = sum(wm_outcome)
	gen prob_uncensored = (N_group-N_censored)/N_group
	gen censnum = prob_uncensored
	sort arm_value pseudo_key tstart
	by arm_value pseudo_key : replace censnum=censnum*censnum[_n-1] if _n!=1
	* generate IPCW
	gen _ipcw = 1 / censdenom
	*
	compress
	save "Sensitivity/SEN5/weight/emulated ipcw `event' SEN5 bs`j'.dta", replace
}
}
*
* Generate KM estimate
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "KM ipcw" _col(80) "SEN5"
	}
	use "Sensitivity/SEN5/weight/emulated ipcw `event' SEN5 bs`j'.dta", clear
	count if invalid == 1
	replace outcome = . if fup_obs != tstop 
	if r(N) == 0 {
		stset tstop [pweight = _ipcw], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(fup bs invalid)
	save "Sensitivity/SEN5/KM ipcw/KM ipcw `event' SEN5 bs`j'.dta", replace
}
}
*
* Finalize bootstrap datasets
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	clear
forvalues j = 0/500 {
	append using "Sensitivity/SEN5/KM ipcw/KM ipcw `event' SEN5 bs`j'.dta"
}
	gen hazard_s_w = 1 - KM_s_w
	gen hazard_ns_w = 1 - KM_ns_w
	gen odds_s_w = hazard_s_w/(1-hazard_s_w)
	gen odds_ns_w = hazard_ns_w/(1-hazard_ns_w)
	gen RR_w = hazard_s_w/hazard_ns_w
	gen diff_w = hazard_s_w - hazard_ns_w
	gen OR_w = odds_s_w / odds_ns_w
	compress
	save "Sensitivity/SEN5/KM ipcw `event' SEN5 bs_all.dta", replace
}
*
cls
* Cumulative incidence (N, %)
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues k = 1/1 {
	use "Sensitivity/SEN5/KM ipcw `event' SEN5 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	keep if _n <= 1 + 500
	
	sum hazard_s_w if bs == 0, d
	scalar hazard_s_mean = r(mean)
	centile hazard_s_w if bs > 0, centile(2.5 97.5)
	scalar hazard_s_cil = r(c_1)
	scalar hazard_s_ciu = r(c_2)
	
	sum hazard_ns_w if bs == 0, d
	scalar hazard_ns_mean = r(mean)
	centile hazard_ns_w if bs > 0, centile(2.5 97.5)
	scalar hazard_ns_cil = r(c_1)
	scalar hazard_ns_ciu = r(c_2)
	
	noi di substr("`event'",1,13) _col(15) "SEN5" _col(30) hazard_s_mean _col(45) hazard_ns_mean
}
}
*
* Absolute risk reduction
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	use "Sensitivity/SEN5/KM ipcw `event' SEN5 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	keep if _n <= 1 + 500
	
	sum diff_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile diff_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "SEN5" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
*
* Relative risk
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	use "Sensitivity/SEN5/KM ipcw `event' SEN5 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	keep if _n <= 1 + 500
	
	sum RR_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile RR_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "SEN5" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
*

******************************************
*** Sen6 - Using IPCW truncated at 1st and 99th percentiles

***Creates directory
capture mkdir "Sensitivity/SEN6"
capture mkdir "Sensitivity/SEN6/KM ipcw_trim"
*
* Generate KM estimate (truncated IPCW)
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
forvalues j = 0/500 {
	if `j'/100 == int(`j'/100) {
	noi di "`event'" _col(25) "subgp_1" _col(40) "bs`j'" _col(55) "KM ipcw_trim" _col(80) "SEN6"
	}
	use "matched_10/subgp1/weight/weight `event' subgp_1 bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = _ipcw_trim], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(fup bs invalid)
	save "Sensitivity/SEN6/KM ipcw_trim/KM `event' SEN6 bs`j'.dta", replace
}
}
*
* Finalize bootstrap datasets
foreach event in covid_hosp MMMI c_section preterm_birth_comb {
clear
forvalues j = 0/500 {
	append using "Sensitivity/SEN6/KM ipcw_trim/KM `event' SEN6 bs`j'.dta"
}
	gen hazard_s_w = 1 - KM_s_w
	gen hazard_ns_w = 1 - KM_ns_w
	gen odds_s_w = hazard_s_w/(1-hazard_s_w)
	gen odds_ns_w = hazard_ns_w/(1-hazard_ns_w)
	gen RR_w = hazard_s_w/hazard_ns_w
	gen diff_w = hazard_s_w - hazard_ns_w
	gen OR_w = odds_s_w / odds_ns_w
	compress
	save "Sensitivity/SEN6/KM ipcw_trim `event' SEN6 bs_all.dta", replace
}
*
cls
* Cumulative incidence (N, %)
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	use "Sensitivity/SEN6/KM ipcw_trim `event' SEN6 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	keep if _n <= 1 + 500
	
	sum hazard_s_w if bs == 0, d
	scalar hazard_s_mean = r(mean)
	centile hazard_s_w if bs > 0, centile(2.5 97.5)
	scalar hazard_s_cil = r(c_1)
	scalar hazard_s_ciu = r(c_2)
	
	sum hazard_ns_w if bs == 0, d
	scalar hazard_ns_mean = r(mean)
	centile hazard_ns_w if bs > 0, centile(2.5 97.5)
	scalar hazard_ns_cil = r(c_1)
	scalar hazard_ns_ciu = r(c_2)
	
	noi di substr("`event'",1,13) _col(15) "SEN6" _col(30) hazard_s_mean _col(45) hazard_ns_mean
}
*
* Absolute risk reduction
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	use "Sensitivity/SEN6/KM ipcw_trim `event' SEN6 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	keep if _n <= 1 + 500
	
	sum diff_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile diff_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "SEN6" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
*
* Relative risk
qui foreach event in covid_hosp MMMI c_section preterm_birth_comb {
	use "Sensitivity/SEN6/KM ipcw_trim `event' SEN6 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	keep if _n <= 1 + 500
	
	sum RR_w if bs == 0, d
	scalar bs_mean = r(mean)
	scalar bs_p50 = r(p50)
	
	centile RR_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "SEN6" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
*

