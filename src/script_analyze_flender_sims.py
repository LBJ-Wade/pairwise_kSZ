import os

################################
# Soergel footprint (1400 deg^2)
################################
# Velocity &  kSZ only
os.system('python analyze_flender_sims.py \
		   -fname_cat ../data/kSZ_tables/T_kSZ_Flender_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP.pkl.gz \
		   -N_bootstraps 300\
		   -do_velocity 1\
		   -fname_output ../data/pairwise/pw_T_kSZ_Flender_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP_0.1_z_0.8_1e14_M200_3e14.pkl.gz &\
			')

# kSZ only + photo-z errs
os.system('python analyze_flender_sims.py \
		   -fname_cat ../data/kSZ_tables/T_kSZ_Flender_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP.pkl.gz \
		   -N_bootstraps 300\
		   -do_velocity 0 \
		   -sigma_photoz 0.01 \
		   -fname_output ../data/pairwise/pw_T_kSZ_Flender_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP_0.1_z_0.8_1e14_M200_3e14_sigmaphz0.01.pkl.gz\
			> logs/pw_T_kSZ_Flender_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP_0.1_z_0.8_1e14_M200_3e14_sigmaphz0.01.pkl.log & ')

# kSZ + CMB + noise
os.system('python analyze_flender_sims.py \
		   -fname_cat ../data/kSZ_tables/T_kSZ_Flender_CMB_noise18muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP.pkl.gz \
		   -N_bootstraps 300 \
		   -do_velocity 0 \
		   -fname_output ../data/pairwise/pw_T_kSZ_Flender_CMB_noise18muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP_0.1_z_0.8_1e14_M200_3e14.pkl.gz \
			> logs/pw_T_kSZ_Flender_CMB_noise18muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP_0.1_z_0.8_1e14_M200_3e14.log &')

# kSZ + tSZ + CMB + noise
os.system('python analyze_flender_sims.py \
		   -fname_cat ../data/kSZ_tables/T_kSZ_tSZ_Flender_CMB_noise18muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP.pkl.gz \
		   -N_bootstraps 300\
		   -do_velocity 0 \
		   -fname_output ../data/pairwise/pw_T_kSZ_tSZ_Flender_CMB_noise18muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP_0.1_z_0.8_1e14_M200_3e14.pkl.gz \
		   > logs/pw_T_kSZ_tSZ_Flender_CMB_noise18muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP_0.1_z_0.8_1e14_M200_3e14.log & ')

# kSZ + tSZ + CMB + noise + photo-z errs 
os.system('python analyze_flender_sims.py \
		   -fname_cat ../data/kSZ_tables/T_kSZ_tSZ_Flender_CMB_noise18muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP.pkl.gz \
		   -N_bootstraps 300\
		   -do_velocity 0 \
		   -sigma_photoz 0.01 \
		   -fname_output ../data/pairwise/pw_T_kSZ_tSZ_Flender_CMB_noise18muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP_0.1_z_0.8_1e14_M200_3e14_sigmaphz0.01.pkl.gz\
			> logs/pw_T_kSZ_tSZ_Flender_CMB_noise18muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP_0.1_z_0.8_1e14_M200_3e14_sigmaphz0.01.log &')

################################
# SPTpol footprint (500 deg^2)
################################

# kSZ + tSZ + CMB + noise + photo-z errs 
os.system('python analyze_flender_sims.py \
		   -fname_cat ../data/kSZ_tables/T_kSZ_tSZ_Flender_CMB_noise5muK_fwhm1.2arcmin_footprintSPTpol_theta2.0arcmin_AP.pkl.gz \
		   -N_bootstraps 300\
		   -do_velocity 0 \
		   -sigma_photoz 0.01 \
		   -fname_output ../data/pairwise/pw_T_kSZ_tSZ_Flender_CMB_noise5muK_fwhm1.2arcmin_footprintSPTpol_theta2.0arcmin_AP_0.1_z_0.8_1e14_M200_3e14_sigmaphz0.01.pkl.gz \
		   > logs/pw_T_kSZ_tSZ_Flender_CMB_noise5muK_fwhm1.2arcmin_footprintSPTpol_theta2.0arcmin_AP_0.1_z_0.8_1e14_M200_3e14_sigmaphz0.01.log &')

################################
# SPT-SZ + SPTpol (1400 deg^2)
################################

# kSZ + tSZ + CMB + noise + photo-z errs 
os.system('python analyze_flender_sims.py \
		   -fname_cat ../data/kSZ_tables/T_kSZ_tSZ_Flender_CMB_noise18muK_noise5muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_footprint2SPTpol_theta2.0arcmin_AP.pkl.gz \
		   -N_bootstraps 300\
		   -do_velocity 0 \
		   -sigma_photoz 0.01 \
		   -fname_output ../data/pairwise/pw_T_kSZ_tSZ_Flender_CMB_noise18muK_noise5muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_footprint2SPTpol_theta2.0arcmin_AP_0.1_z_0.8_1e14_M200_3e14_sigmaphz0.01.pkl.gz \
		   > logs/pw_T_kSZ_tSZ_Flender_CMB_noise18muK_noise5muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_footprint2SPTpol_theta2.0arcmin_AP_0.1_z_0.8_1e14_M200_3e14_sigmaphz0.01.log &')

################################
# SPT-3G        (1400 deg^2)
################################

# kSZ + tSZ + CMB + noise + photo-z errs 
os.system('python analyze_flender_sims.py \
		   -fname_cat ../data/kSZ_tables/T_kSZ_tSZ_Flender_CMB_noise5muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP.pkl.gz \
		   -N_bootstraps 300\
		   -do_velocity 0 \
		   -sigma_photoz 0.01 \
		   -fname_output ../data/pairwise/pw_T_kSZ_tSZ_Flender_CMB_noise5muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP_0.1_z_0.8_1e14_M200_3e14_sigmaphz0.01.pkl.gz\
			> logs/pw_T_kSZ_tSZ_Flender_CMB_noise5muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP_0.1_z_0.8_1e14_M200_3e14_sigmaphz0.01.log &')



