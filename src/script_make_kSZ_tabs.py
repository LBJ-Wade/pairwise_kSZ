import os

os.system('python make_T_kSZ_table_flender_sims.py \
				-AP_type AVG \
				-theta_arcmin 2.0\
				-beam 1.2 \
				-footprint SPT-SZ \
				-fname_output ../data/kSZ_tables/T_kSZ_Flender_fwhm1.2arcmin_footprintSPT-SZ_theta2.0arcmin_AVG.pkl.gz' )

os.system('python make_T_kSZ_table_flender_sims.py \
				-AP_type AP \
				-theta_arcmin 2.0\
				-beam 1.2 \
				-footprint SPT-SZ \
				-fname_output ../data/kSZ_tables/T_kSZ_Flender_fwhm1.2arcmin_footprintSPT-SZ_theta2.0arcmin_AP.pkl.gz' )

### 

os.system('python make_T_kSZ_table_flender_sims.py \
				-AP_type AP \
				-theta_arcmin 2.0\
				-beam 1.2 \
				-footprint SPTpol \
				-fname_output ../data/kSZ_tables/T_kSZ_Flender_fwhm1.2arcmin_footprintSPTpol_theta2.0arcmin_AP.pkl.gz' )


os.system('python make_T_kSZ_table_flender_sims.py \
				-AP_type AP \
				-theta_arcmin 2.0\
				-footprint SPT-SZ \
				-beam 1.2 \
				-tSZ True \
				-fname_output ../data/kSZ_tables/T_kSZ_tSZ_Flender_fwhm1.2arcmin_footprintSPT-SZ_theta2.0arcmin_AP.pkl.gz' )

os.system('python make_T_kSZ_table_flender_sims.py \
				-AP_type AP \
				-theta_arcmin 2.0\
				-footprint SPT-SZ \
				-beam 1.2 \
				-tSZ True \
				-cmb True \
				-fname_output ../data/kSZ_tables/T_kSZ_tSZ_Flender_CMB_fwhm1.2arcmin_footprintSPT-SZ_theta2.0arcmin_AP.pkl.gz' )

os.system('python make_T_kSZ_table_flender_sims.py \
				-AP_type AP \
				-theta_arcmin 2.0\
				-footprint SPT-SZ \
				-beam 1.2 \
				-cmb True \
				-noise 18. \
				-fname_output ../data/kSZ_tables/T_kSZ_Flender_CMB_noise18muK_fwhm1.2arcmin_footprintSPT-SZ_theta2.0arcmin_AP.pkl.gz' )

os.system('python make_T_kSZ_table_flender_sims.py \
				-AP_type AP \
				-theta_arcmin 2.0\
				-footprint SPT-SZ \
				-beam 1.2 \
				-tSZ True \
				-cmb True \
				-noise 18. \
				-fname_output ../data/kSZ_tables/T_kSZ_tSZ_Flender_CMB_noise18muK_fwhm1.2arcmin_footprintSPT-SZ_theta2.0arcmin_AP.pkl.gz' )

os.system('python make_T_kSZ_table_flender_sims.py \
				-AP_type AP \
				-theta_arcmin 2.0\
				-footprint SPT-SZ \
				-footprint2 SPTpol \
				-beam 1.2 \
				-tSZ True \
				-cmb True \
				-noise 18. \
				-noise2 5. \
				-fname_output ../data/kSZ_tables/T_kSZ_tSZ_Flender_CMB_noise18muK_noise5muK_fwhm1.2arcmin_footprintSPT-SZ_footprint2SPTpol_theta2.0arcmin_AP.pkl.gz' )
########
# Soergel
########

# kSZ-only
os.system('python make_T_kSZ_table_flender_sims.py \
				-AP_type AP \
				-theta_arcmin 2.0\
				-beam 1.2 \
				-footprint SPT-SZ_Soergel \
				-fname_output ../data/kSZ_tables/T_kSZ_Flender_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP.pkl.gz' )

# kSZ + tSZ + CMB + noise
os.system('python make_T_kSZ_table_flender_sims.py \
				-AP_type AP \
				-theta_arcmin 2.0\
				-footprint SPT-SZ_Soergel \
				-beam 1.2 \
				-tSZ True \
				-cmb True \
				-noise 18. \
				-fname_output ../data/kSZ_tables/T_kSZ_tSZ_Flender_CMB_noise18muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP.pkl.gz' )


# kSZ + tSZ + CMB + noise (SPT-3G)
os.system('python make_T_kSZ_table_flender_sims.py \
				-AP_type AP \
				-theta_arcmin 2.0\
				-footprint SPT-SZ_Soergel \
				-beam 1.2 \
				-tSZ True \
				-cmb True \
				-noise 5. \
				-fname_output ../data/kSZ_tables/T_kSZ_tSZ_Flender_CMB_noise5muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP.pkl.gz' )



os.system('python make_T_kSZ_table_flender_sims.py \
				-AP_type AP \
				-theta_arcmin 2.0\
				-footprint SPT-SZ_Soergel \
				-beam 1.2 \
				-tSZ True \
				-cmb True \
				-noise 18. \
				-fname_output ../data/kSZ_tables/T_kSZ_tSZ_Flender_CMB_noise18muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP.pkl.gz' )

# kSZ + CMB + noise
os.system('python make_T_kSZ_table_flender_sims.py \
				-AP_type AP \
				-theta_arcmin 2.0\
				-footprint SPT-SZ_Soergel \
				-beam 1.2 \
				-cmb True \
				-noise 18. \
				-fname_output ../data/kSZ_tables/T_kSZ_Flender_CMB_noise18muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_theta2.0arcmin_AP.pkl.gz' )

# kSZ + tSZ + CMB + noise (SPT-SZ + SPTpol)
os.system('python make_T_kSZ_table_flender_sims.py \
				-AP_type AP \
				-theta_arcmin 2.0\
				-footprint SPT-SZ_Soergel \
				-footprint2 SPTpol \
				-beam 1.2 \
				-tSZ True \
				-cmb True \
				-noise 18. \
				-noise2 5. \
				-fname_output ../data/kSZ_tables/T_kSZ_tSZ_Flender_CMB_noise18muK_noise5muK_fwhm1.2arcmin_footprintSPT-SZ_Soergel_footprint2SPTpol_theta2.0arcmin_AP.pkl.gz' )

#########
# SPTpol
#########

os.system('python make_T_kSZ_table_flender_sims.py \
				-AP_type AP \
				-theta_arcmin 2.0\
				-footprint SPTpol \
				-beam 1.2 \
				-tSZ True \
				-cmb True \
				-noise 5. \
				-fname_output ../data/kSZ_tables/T_kSZ_tSZ_Flender_CMB_noise5muK_fwhm1.2arcmin_footprintSPTpol_theta2.0arcmin_AP.pkl.gz' )

############
# SPTpoldeep
############

os.system('python make_T_kSZ_table_flender_sims.py \
				-AP_type AP \
				-theta_arcmin 2.0\
				-footprint SPTpol_deep \
				-beam 1.2 \
				-tSZ True \
				-cmb True \
				-noise 5. \
				-fname_output ../data/kSZ_tables/T_kSZ_tSZ_Flender_CMB_noise5muK_fwhm1.2arcmin_footprintSPTpol_deep_theta2.0arcmin_AP.pkl.gz' )
