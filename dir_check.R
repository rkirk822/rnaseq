
# dir_check.R
# 
# Make sure a given directory exists, and if there's no slash at the end of the path, add it.

dir_check = function(dirname) {

    # Check that directory exists
    if ( !dir.exists(dirname) ) {stop(paste('Directory', dirname, 'does not exist.'))}
    
    # Make sure it has a / at the end
    if ( ! substr(dirname, nchar(dirname), nchar(dirname))=='/') { dirname = paste(dirname, '/', sep='') }
    
    return(dirname)
    
}
