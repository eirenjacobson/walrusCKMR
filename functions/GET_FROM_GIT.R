"GET_FROM_GIT" <-
function(){
  dirz <- cq( functions, simulation, results, scripts)
  unlink( dirz)
  for( adir in dirz){
    file.copy( file.path( 'd:/github/eirenj/walrusCKMR', adir), 
        adir, recursive=TRUE)
  }
  
  # Load the functions (also in main script)
  lapply( list.files( './functions', full=TRUE), source) 
}
