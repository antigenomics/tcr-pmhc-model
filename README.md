# tcr-pmhc-model
Scripts for modelling TCR:pMHC structures

``ch_scripts/ch_select_models/`` - initialization/preprocessing scripts

* ``ch_remove_chains.py`` - remove certain chains (TRA/TRB/antigen/MHC) from structures
* ``ch_analyse_vdjdb.py`` - select VDJdb records that have matching (by V segment, J segment and CDR3 length) templates in structural data

``modeller/ch_run_modeller.py`` - execute modeller
