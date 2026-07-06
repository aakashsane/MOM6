This is my MOM6 that has the new module to prescribe flux of tracer at the surface in a lat-lon rectangle patch.


The MOM6-examples git verison is:

<!--- 
[as8486@stellar-amd MOM6-examples]$ git branch
* dev/gfdl
[as8486@stellar-amd MOM6-examples]$ git log
commit 5abb55d09e31ef95e4e35299a1f09678798d3214 (HEAD -> dev/gfdl, origin/gfdl-to-main-2025-09-25, origin/dev/gfdl, origin/HEAD)
Author: Marshall Ward <marshall.ward@noaa.gov>
Date:   Thu Sep 25 14:30:15 2025 -0400

    SIS2: Neglect massless defaults (#223)
    
    - NOAA-GFDL/SIS2@4cbe47c Neglect massless defaults (#223)
    - NOAA-GFDL/SIS2@e9fa2e9  Adding brine plume interface (#222)
    - NOAA-GFDL/SIS2@c4f3406 Correct overly long lines and trailing whitespace
    - NOAA-GFDL/SIS2@378c777 +Add set_segnum_signs and OBC loop ranges on PE
    - NOAA-GFDL/SIS2@b17231f +Signed OBC%segnum_u and OBC%segnum_v in SIS2
    - NOAA-GFDL/SIS2@4f86208 Fixed a bug that could incorrectly trigger a warning about non-conservation of ice thickness and concentration during the distribution routine

commit ba6a3a644caf4ec9daf7890e194f207875374b5c
Author: MOM6 bot <mom6bot@users.noreply.github.com>
Date:   Wed Sep 24 14:11:51 2025 -0400

    MOM6: +Add the new parameter RESOLN_FUNCTION_OBC_BUG
    
    - NOAA-GFDL/MOM6@39ab7d508 +Add the new parameter RESOLN_FUNCTION_OBC_BUG
    - NOAA-GFDL/MOM6@665c760d8 Flip the order of acceleration and velocity chksum

--->

