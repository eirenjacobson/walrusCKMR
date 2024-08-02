"make_n_comps" <- function( 
    m_YA,
    m_F_YANL,
    m_F_YAself,
    m_F_YD,
    FIRST_PDYEAR,
    AMAT,
    MAX_BGAP_HSPS
){
  ## Abstracted from add_data, for re-use in pure designs
  SYEARS <- dimseq( m_YA)[[1]]
  AGES <- dimseq( m_YA)[[2]]
  NURSINGNESS <- dimseq( m_F_YANL)[[3]]
  LETHALITY <- dimseq(m_F_YANL)[[4]] 
  
  SELF1AGES <- dimseq( m_F_YAself)[[2]]
  DEVSTAGES <- dimseq( m_F_YD)[[2]]
  
  # _Could_ do consistency checks on all input dimensions
  # but leave that to R (it's really up to the caller!!!)
  
  # Cheekily, I am really using A (exact age) but calling it E
  # EKJ 2024-06-21 changed this to A to avoid confusion
  n_comp_MOP_AYNL <- autoloop(
    aj= AGES, yj = SYEARS,
    ac = AGES, yc = SYEARS, nc=NURSINGNESS, lc = LETHALITY, {
      bad <- yc - ac
      bju <- yj - aj
      
        (bju >= FIRST_PDYEAR) *
        (bju >= (bad+AMAT)) *
        (yj != yc) * # should only apply to pre-max-weanage, but...
        m_YA[yj, aj] * m_F_YANL[yc, ac, nc, lc] # m_F_YANL
    })
  n_comp_MOP_AYL <- sumover( n_comp_MOP_AYNL, 'nc')
  

  n_comp_XmHSP_AY <- autoloop(
    a1= AGES, y1 = SYEARS,
    a2 = AGES, y2 = SYEARS, {
      b1 <- y1 - a1
      b2 <- y2 - a2
      
        (b1 >= FIRST_PDYEAR) * (b2 >= FIRST_PDYEAR) *
        (b2>b1) * # don't double-count
        ((b2-b1) <= MAX_BGAP_HSPS) * # GGP risk beyond this
        m_YA[y1,a1] * m_YA[y2,a2]
    })
  
  n_comp_selfP_YADY <- autoloop(
    y1= SYEARS, a1= SELF1AGES, # add lethality here?
    d2= DEVSTAGES, y2= SYEARS,{
      (y2 > y1) * # duhhh
        # keep m_FLive_YAself, option to add lethality to m_F_YD
        m_F_YAself[ y1, a1] * m_F_YD[ y2, d2]
    })
  
  # Return all n_comp_<blah>
  returnList( mget( ls( environment(), pattern='^n_comp_')))
}
