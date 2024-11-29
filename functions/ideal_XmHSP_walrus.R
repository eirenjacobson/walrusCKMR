"ideal_XmHSP_walrus" <-
function( nlocal=sys.parent()) mlocal({
## Non-spatial only here.
# Spatial version has the extra complexity that juve will stay with mum
# until weaning (by defn) so the two sibs are not moving independently
  
  Pr_XmHSP_bb <- autoloop(
      b1= PDYEARS, b2= PDYEARS,
    (b1 < b2) *   # non-symmn
      Pr_fadsurv_t[ abs( b2-b1)] * # avoid OOB
      Pr_breedyagain_Db[ abs( b2-b1)] * 
      recip_Nfad_y[ pmax( b1, b2)] *
      recip_ppn_breedy
    )
})
<bytecode: 0x000002212e385488>
