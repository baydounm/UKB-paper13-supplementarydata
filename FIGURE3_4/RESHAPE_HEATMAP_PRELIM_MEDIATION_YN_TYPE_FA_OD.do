**Mediation results with 144 proteins: reshaping of file**
**Copy paste from excel, each of FA, MD, ICVF, ISOVF, and OD results*
**Run the following reshape command and then send to Yi-Han to do the heat map**



cd "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER8F_LE8PROT_DMRI\MANUSCRIPT_FINAL\GITHUB\FIGURES\FIGURE3_4"

**********************************************FA*****************************************************
**NO MEDIATION FINDGINS**
use NO_MEDIATION_GROUPA_FA,clear


capture replace p="0" if p=="<0.001"
destring p, replace

save NO_MEDIATION_GROUPA_FAfin,replace


reshape wide beta se z p lcl ucl, i(protein) j(fourwaydecomp, string)

save NO_MEDIATION_GROUPA_FAwide, replace




**INCONSISTENT MEDIATION FINDINGS**
use INCONSISTENT_MEDIATION_GROUPB_FA,clear


capture replace p="0" if p=="<0.001"
destring p, replace

save INCONSISTENT_MEDIATION_GROUPB_FAfin,replace


reshape wide beta se z p lcl ucl, i(protein) j(fourwaydecomp, string)

save INCONSISTENT_MEDIATION_GROUPB_FAwide, replace



**CONSISTENT MEDIATION FINDINGS**
use CONSISTENT_MEDIATION_GROUPC_FA,clear


capture replace p="0" if p=="<0.001"
destring p, replace

save CONSISTENT_MEDIATION_GROUPC_FAfin,replace


reshape wide beta se z p lcl ucl, i(protein) j(fourwaydecomp, string)

save CONSISTENT_MEDIATION_GROUPC_FAwide, replace



***********************************************OD****************************************************

**NO OR INCONSISTENT MEDIATION FINDINGS**" 
use INCONSISTENT_MEDIATION_GROUPAB_OD,clear


capture replace p="0" if p=="<0.001"
destring p, replace

save INCONSISTENT_MEDIATION_GROUPAB_ODfin,replace


reshape wide beta se z p lcl ucl, i(protein) j(fourwaydecomp, string)

save INCONSISTENT_MEDIATION_GROUPAB_ODwide, replace



**CONSISTENT MEDIATION FINDINGS**
use CONSISTENT_MEDIATION_GROUPC_OD,clear

replace p="0" if p=="<0.001"
destring p, replace

save CONSISTENT_MEDIATION_GROUPC_ODfin, replace


reshape wide beta se z p lcl ucl, i(protein) j(fourwaydecomp, string)

save CONSISTENT_MEDIATION_GROUPC_ODwide, replace



