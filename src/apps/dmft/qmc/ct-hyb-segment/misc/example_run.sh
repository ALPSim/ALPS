module load alps_cthyb

npts_tau=256
beta=5.0

python gen_sample_hybridization.py --gamma 0.2 --npts ${npts_tau} --beta ${beta}
alps_cthyb --N_ORBITALS 2 --U 4 --MU 2 --N_TAU ${npts_tau} --BETA ${beta} --DELTA delta_tau.dat --N_MEAS 3000 --THERMALIZATION 100 --SWEEPS 200000 --TEXT_OUTPUT 1
