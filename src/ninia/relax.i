
f"""
 &CONTROL
    calculation = 'relax',
    prefix      = '{prefix}',
    outdir      = '{output_dir}',
    pseudo_dir  = '{pseudo_dir}',
 /
 &SYSTEM
    ibrav     = 0,
    nat  = {num_atoms},
    ntyp = {num_elem},
    ecutwfc = {ecutwfc},
    input_dft = '{functional}',

    occupations = 'smearing',
    smearing = 'mv',
    degauss = 0.01,
 /
 &ELECTRONS
    conv_thr = {conv_thr} # 1.0d-8
    mixing_beta = {mixing_beta} # 0.15d0
 /
 &IONS
 /
 &CELL
 /

ATOMIC_SPECIES
{atomic_species}

#   Ru  101.07   Ru_ONCV_PBE-1.0.oncvpsp.upf
#   C   12.011   C.pbe-n-kjpaw_psl.1.0.0.UPF
#   H   1.008    H_ONCV_PBE-1.0.oncvpsp.upf

CELL_PARAMETERS angstrom
{cell_parameters}

#9.36443933730765   5.40656161739439   0.00000000000000
#0.00000000000000   5.40656161739439   0.00000000000000
#0.00000000000000   0.00000000000000   18.621657682203

ATOMIC_POSITIONS angstrom
{atomic_positions}

K_POINTS automatic
{k_points}

#   4 4 4 0 0 0


"""