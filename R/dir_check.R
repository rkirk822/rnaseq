#' Check directory existence
#'
#' Make sure a given directory exists, and if there's no slash at the end of the path, add it.
#' @param dirname String - directory to check
#' @return String - directory with slash at end of name (if it didn't have one already)
#' @details Causes error if specified directory doesn't exist.  Puts slash at end of name if it's not there already.
#' @author Emma Myers
#' @export

dir_check = function(dirname) {

    # Check that directory exists
    if ( !dir.exists(dirname) ) {stop(paste('Directory', dirname, 'does not exist.'))}

    # Make sure it has a / at the end
    if ( ! substr(dirname, nchar(dirname), nchar(dirname))=='/') { dirname = paste(dirname, '/', sep='') }

    return(dirname)

}
