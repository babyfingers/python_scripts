$molecule
xxCHARGExx xxMULTxx
xxGEOMETRYxx
$end

$rem
jobtype             	sp
exchange	    	omegaB97X
basis           	6-31G*
mem_static 		240
mem_total 		2000
symmetry        	false
sym_ignore 		true
cis_n_roots		6
cis_singlets		true
cis_triplets		false
scf_guess               read
cis_guess_disk          true
cis_guess_disk_type     2
xxREM_LD_TYPExx      	xxLD_NSTATExx
xxREM_COMMANDxx$end

$localized_diabatization
xxLD_COMMENTxx
xxLD_STATESxx
$end



