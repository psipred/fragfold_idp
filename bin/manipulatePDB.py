import numpy as np
import random
import itertools
import glob
import re

"""
    manipulate PDB files
        * read PDB
        * save PDB

        * generate ensemble statistics
        * calculate ensemble statistics
"""


def ReadPDB(filename, chain_mode='CA'):
    '''
    reads a PDB file with single structure or an ensemble of structures
    (either main-chain or CA-only or main-chain+CB)

    input:
    filename incl. path,
    *optional: specify if main-chain [MC] is to be read,
    main-chain+CB [CB] or [CA] only [default]

    returns:
    numpy.array (or arrays) of residues and conformations
    the order of output N, CA, CB, C, O
    order of info in array - [conformation,amino-acid,xyz]
    '''

# pre-processing

    if chain_mode != 'MC' and chain_mode != 'CA' and chain_mode != 'CB':
        raise ValueError('Chain mode in ReadPDB() incorrectly specified! \
                         Specify either CA, CB or MC')

# read file
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

# setup variables
    conf = 0
    aa = 0
    max_aa = 0
    at_count = 0

# pre-allocate numpy array/s
# max. 500 residues in *MAX_N_CONF* conformations
    max_n_conf = 500
    CA_xyz = np.zeros((max_n_conf, 500, 3))
    if chain_mode == 'MC' or chain_mode == 'CB':
        C_xyz = np.zeros((max_n_conf, 500, 3))
        N_xyz = np.zeros((max_n_conf, 500, 3))
        O_xyz = np.zeros((max_n_conf, 500, 3))
    if chain_mode == 'CB':
        CB_xyz = np.zeros((max_n_conf, 500, 3))

# read PDB into numpy array
    for line in lines:
        if line.startswith('ATOM') and 'CA' in line:
            CA_xyz[conf, aa, :] = re.sub('-', ' -', line[29:54]).split()
            at_count += 1
        elif line.startswith('ATOM') and ' C ' in line[0:20] and \
                            ('MC' in chain_mode or 'CB' in chain_mode):
            C_xyz[conf, aa, :] = re.sub('-', ' -', line[29:54]).split()
            at_count += 1
        elif line.startswith('ATOM') and ' N ' in line[0:20] and \
                            ('MC' in chain_mode or 'CB' in chain_mode):
            N_xyz[conf, aa, :] = re.sub('-', ' -', line[29:54]).split()
            at_count += 1
        elif line.startswith('ATOM') and ' O ' in line[0:20] and \
                            ('MC' in chain_mode or 'CB' in chain_mode):
            O_xyz[conf, aa, :] = re.sub('-', ' -', line[29:54]).split()
            at_count += 1
        elif line.startswith('ATOM') and 'CB' in chain_mode and 'CB' in line:
            CB_xyz[conf, aa, :] = re.sub('-', ' -', line[29:54]).split()
            at_count += 1
        if at_count == 4 and 'GLY' in line and 'CB' in chain_mode:
            CB_xyz[conf, aa, :] = GenerateCB(CA_xyz[conf, aa, :],
                                             N_xyz[conf, aa, :],
                                             C_xyz[conf, aa, :])
            at_count += 1

        if at_count == 1 and 'CA' in chain_mode:
            aa += 1
            at_count = 0
        elif at_count == 4 and 'MC' in chain_mode:
            aa += 1
            at_count = 0
        elif at_count == 5 and 'CB' in chain_mode:
            aa += 1
            at_count = 0

        if line.startswith('END\n') or \
           line.startswith('TER') \
           or line.startswith('END '):
            if aa > 0:
                conf += 1
            if aa > max_aa:
                max_aa = aa
            aa = 0

# produce output
    CA_xyz = CA_xyz[:conf, :max_aa, :]
    if chain_mode == 'MC' or chain_mode == 'CB':
        C_xyz = C_xyz[:conf, :max_aa, :]
        N_xyz = N_xyz[:conf, :max_aa, :]
        O_xyz = O_xyz[:conf, :max_aa, :]
    if chain_mode == 'CB':
        CB_xyz = CB_xyz[:conf, :max_aa, :]

# return
    if conf == 1:
        CA_xyz = CA_xyz[0]
        if chain_mode == 'CB' or chain_mode == 'MC':
            C_xyz = C_xyz[0]
            N_xyz = N_xyz[0]
            O_xyz = O_xyz[0]
        if chain_mode == 'CB': CB_xyz = CB_xyz[0]
    if chain_mode == 'CA':
        return CA_xyz
    elif chain_mode == 'MC':
        return N_xyz, CA_xyz, C_xyz, O_xyz
    elif chain_mode == 'CB':
        return N_xyz, CA_xyz, CB_xyz, C_xyz, O_xyz


def SavePDB(xyz_CA, residues, target):

    import datetime.date
    formatPDB = "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n"

    temp_f = 10.0
    x_ca = []
    y_ca = []
    z_ca = []
#   x_cb = []
#   y_cb = []
#   z_cb = []

    for coords_CA in xyz_CA:  # split xyzCA!!!
        x_ca.append(coords_CA.split()[0])
        y_ca.append(coords_CA.split()[1])
        z_ca.append(coords_CA.split()[2])
#        x_cb.append(coords_CB[0])
#        y_cb.append(coords_CB[1])
#        z_cb.append(coords_CB[2])

    res_name = []
    for res in residues:
        res_name.append(res.replace('\n', ''))

    f = open(target, 'w')

    at_no = 1
    res_no = 0
    at_name = ' CA '  # !!!!!!!!

    if header_info:
        f.write('HEADER    %s') % header_info

    f.write('EXPDTA    GENERATED USING PYTHON BY TOMASZ KOSCIOLEK ON %s') %
            datetime.date.today().strftime('%d/%m/%Y')

    for line in zip(res_name, x_ca, y_ca, z_ca):
        res_no += 1

        f.write(str(formatPDB) % (at_no, at_name, line[0], res_no,
                float(line[1]),
                float(line[2]),
                float(line[3]),
                temp_f))
        at_no += 1

    f.write('END')
    f.close()


def GenerateEnsStats(nmr_inp, ff_inp, FF_ens_size, FF_n_ens, SW_size, name_pickle):

	##########################
	##	input files	##
	##    1) nmr_inp	##
	##    2) ff_inp		##
	##########################

	#############
	# NMR input #
	#############

	## READ PDB FILES
	NMR_PDB_ens = ReadPDB(nmr_inp) # mode CA

	## FIND NUMBER OF CONFORMATIONS IN THE FILE
	NMR_n_conf = NMR_PDB_ens.shape[0]
	NMR_n_ens = 1
	prot_len = NMR_PDB_ens.shape[1]

	## GENERATE *PER-RESIDUE SLIDING WINDOW RMSD*
	## STORE ARRAY (VARIABLE)
	NMR_PerRes_RMSD = Sliding_Window(NMR_PDB_ens, sw = int(SW_size))

	## INTER-RESIDUE RMSD FULL ALIGNMENT - **NO** SLIDING WINDOW!
	NMR_InterRes_RMSD = np.zeros((NMR_n_ens, prot_len))
	for ens in range(NMR_n_ens):
		NMR_InterRes_RMSD[ens,:] = multi_RMSD(NMR_PDB_ens, per_res = True, max_iter = 6, thr_val = 0.2)

	## CALCULATE RMSD FOR EACH PAIR OF MODELS IN AN ENSEMBLE
	NMR_tot_RMSD = np.zeros((NMR_n_ens, NMR_n_conf*(NMR_n_conf-1)/2))
	for ens in range(NMR_n_ens):
		model_iter = 0
		for confA in range(NMR_n_conf):
		 for confB in range(NMR_n_conf):
		  if confA > confB:
		   NMR_tot_RMSD[ens,model_iter] = multi_RMSD(NMR_PDB_ens[confA,:,:], NMR_PDB_ens[confB,:,:], per_res = False, max_iter = 3, thr_val = 1.0)
		   model_iter += 1

	#############
	# FF  input #
	#############

	## READ *ENS* FILES
	if '*' in ff_inp:
		FF_ens = []
		for inp_file in glob.glob(ff_inp):
			FF_ens.append(ReadPDB(inp_file))
		FF_ens = np.array(FF_ens)
	else:
		FF_ens = ReadPDB(ff_inp) # in CA mode

	## GENERATE 100 RANDOM ENSEMBLES OF 10 STRUCTURES
	n_conf = FF_ens.shape[0]
	prot_len = FF_ens.shape[1]

	conf_list = [i for i in range(n_conf)]
	FF_ens_list = [random.sample(conf_list, FF_ens_size) for i in xrange(FF_n_ens)]

	## CALCULATE *PER-RESIDUE SLIDING-WINDOW RMSD* FOR EACH ENSEMBLE IN FF_ens_list
	FF_PerRes_RMSD = np.zeros((FF_n_ens, prot_len))			# np.array
	for ens in enumerate(FF_ens_list):
		FF_PerRes_RMSD[ens[0],:] = Sliding_Window(FF_ens[ens[1]], sw = int(SW_size))

	## INTER-RESIDUE RMSD FULL ALIGNMENT - **NO** SLIDING WINDOW!
	FF_InterRes_RMSD = np.zeros((FF_n_ens, prot_len))			# np.array
	for ens in enumerate(FF_ens_list):
		FF_InterRes_RMSD[ens[0],:] = multi_RMSD(FF_ens[ens[1]], per_res = True, max_iter = 6, thr_val = 0.2)


	## TOTAL RMSD FOR EACH PAIR OF MODELS IN AN ENSEMBLE
	FF_tot_RMSD = np.zeros((FF_n_ens, FF_ens_size*(FF_ens_size-1)/2))	# np.array
	for ens in range(FF_n_ens):
		model_iterator = 0
		for FF_conf_A, FF_conf_B in itertools.product(FF_ens_list[ens],FF_ens_list[ens]):
			if FF_conf_A != FF_conf_B and (FF_conf_A > FF_conf_B):
				FF_tot_RMSD[ens,model_iterator] = multi_RMSD(FF_ens[FF_conf_A], FF_ens[FF_conf_B], per_res=False, max_iter=3, thr_val = 1.0)
				model_iterator += 1

	## FF ENERGY FOR EACH MODEL IN AN ENSEMBLE
	FF_tot_E = np.zeros((FF_n_ens, FF_ens_size))				# np.array
	f = open(ff_inp, 'r')
	lines = f.readlines()
	f.close()
	FF_E_list = np.zeros(n_conf)
	model_iterator = 0
	for line in lines:
		if line.startswith('HEADER') and 'Potential' in line:
			FF_E_list[model_iterator] = float(line.split()[3])
			model_iterator += 1
	for ens in range(FF_n_ens):
		model_iterator = 0
		for model_id in FF_ens_list[ens]:
			FF_tot_E[ens, model_iterator] = FF_E_list[model_id]
			model_iterator += 1

	if name_pickle:
		# n_n, n_tot, n_pr, n_ir, ff_n, ff_list, ff_tot, ff_pr, ff_ir, ff_e
		DumpVars(name_pickle, NMR_n_ens, NMR_tot_RMSD, NMR_PerRes_RMSD, NMR_InterRes_RMSD, FF_n_ens, FF_ens_list, FF_tot_RMSD, FF_PerRes_RMSD, FF_InterRes_RMSD, FF_tot_E)

	return CalculateEnsStats(NMR_n_ens, NMR_tot_RMSD, NMR_PerRes_RMSD, NMR_InterRes_RMSD, FF_n_ens, FF_tot_RMSD, FF_PerRes_RMSD, FF_InterRes_RMSD, FF_tot_E)

def CalculateEnsStats(NMR_n_ens, NMR_tot_RMSD, NMR_PerRes_RMSD, NMR_InterRes_RMSD, FF_n_ens, FF_tot_RMSD, FF_PerRes_RMSD, FF_InterRes_RMSD, FF_tot_E):
	#################
	#### RESULTS ####
	#################

	### allocate numpy arrays
	# Rs
	res_Rs = np.zeros(FF_n_ens)

	# FF RMSD
	res_FF_uRMSD = np.zeros(FF_n_ens)
	res_FF_sdRMSD = np.zeros(FF_n_ens)
	res_FF_uIrRMSD = np.zeros(FF_n_ens)
	res_FF_medRMSD = np.zeros(FF_n_ens)
	res_FF_deltaRMSD = np.zeros(FF_n_ens)

	# FF E
	res_FF_uE = np.zeros(FF_n_ens)
	res_FF_sdE = np.zeros(FF_n_ens)
	res_FF_deltaE = np.zeros(FF_n_ens)
	res_FF_medE = np.zeros(FF_n_ens)

	# NMR RMSD
	res_NMR_uRMSD = np.zeros(NMR_n_ens)
	res_NMR_sdRMSD = np.zeros(NMR_n_ens)
	res_NMR_uIrRMSD = np.zeros(NMR_n_ens)
	res_NMR_medRMSD = np.zeros(NMR_n_ens)
	res_NMR_deltaRMSD = np.zeros(NMR_n_ens)


	### Rs
	for ens in range(FF_n_ens):
		res_Rs[ens] = Rs(NMR_PerRes_RMSD, FF_PerRes_RMSD[ens])

#### RMSD ####

########
## FF ##
########
	### mean RMSD
	for ens in range(FF_n_ens):
		res_FF_uRMSD[ens] = np.mean(FF_tot_RMSD[ens])
	### RMSD SD
		res_FF_sdRMSD[ens] = np.std(FF_tot_RMSD[ens]) 	#### I HAVE TO CONSIDER ONLY ACTUAL *UNIQUE* VALUES
	### RMSD spread (deltaRMSD)
		res_FF_deltaRMSD[ens] = max(FF_tot_RMSD[ens])-min(FF_tot_RMSD[ens])
	### median RMSD
		res_FF_medRMSD[ens] = np.median(FF_tot_RMSD[ens])
	### mean inter-residue RMSD
		res_FF_uIrRMSD[ens] = np.mean(FF_InterRes_RMSD[ens])

#########
## NMR ##
#########
	### mean RMSD
	for ens in range(NMR_n_ens):
		res_NMR_uRMSD[ens] = np.mean(NMR_tot_RMSD[ens])
	### RMSD SD
		res_NMR_sdRMSD[ens] = np.std(NMR_tot_RMSD[ens])
	### mean inter-residue RMSD
		res_NMR_uIrRMSD[ens] = np.mean(NMR_InterRes_RMSD[ens])
	### RMSD spread
		res_NMR_deltaRMSD[ens] = max(NMR_tot_RMSD[ens])-min(NMR_tot_RMSD[ens])
	### median RMSD
		res_NMR_medRMSD[ens] = np.median(NMR_tot_RMSD[ens])

##############

#### Energy ####
	### mean E
	for ens in range(FF_n_ens):
		res_FF_uE[ens] = np.mean(FF_tot_E[ens])
	### E SD
		res_FF_sdE[ens] = np.std(FF_tot_E[ens])
	### energy spread (deltaE)
		res_FF_deltaE[ens] = max(FF_tot_E[ens])-min(FF_tot_E[ens])
	### median E
		res_FF_medE[ens] = np.median(FF_tot_E[ens])
################


### RETURN ###
	return res_Rs, res_FF_uRMSD, res_FF_sdRMSD, res_FF_uIrRMSD, res_FF_medRMSD, res_FF_deltaRMSD, res_FF_uE, res_FF_sdE, res_FF_deltaE, res_FF_medE, res_NMR_uRMSD, res_NMR_sdRMSD, res_NMR_uIrRMSD, res_NMR_medRMSD, res_NMR_deltaRMSD
##############
