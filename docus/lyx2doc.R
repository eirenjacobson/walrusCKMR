# wrapper around lyxport for those not afraid of the keyboard
lick_my_doc <- function(file){
  # first convert to lyxgz ("lick giz"?)
  cmd <- paste("lyx --export lyxgz", file)
  r <- system(cmd)
  
  if(r==0){
    lyxport::lyxzip2word(file)
  }else{
    stop("oh no, something bad happened!")
  }
}

lick_my_doc("docus/paper-ckwaldes.lyx")
