"PUT_TO_GIT" <-
function(){
  # Check if code has be
  funcs <- find.funs( .GlobalEnv)
  sfile <- do.on( funcs, .GlobalEnv[[.]]@srcref@srcfile$filename)
  funcs <- funcs[ sfile=='dummyfile']
  for( f in funcs){
    sink( sprintf( './functions/%s.R', f))
    scatn( '"%s" <-', f)
    print( .GlobalEnv[[f]])
    sink()
  }
  
  file.copy( file.path( 'functions', funcs),
    'd:/github/eirenj/walrusCKMR/functions',
    overwrite=TRUE)
    
  # Scripts:
  dotR <- dir( patt='.[rR]$')
  file.copy( dotR, 'd:/github/eirenj/walrusCKMR', 
    overwrite=TRUE)
}
#<bytecode: 0x0000027765199548>
