#' Check file existence / extension
#'
#' Checks whether a file exists and returns a message that depends on whether that file was supposed to exist.
#' It can also check whether the file has a given extension, and return a message based on that too.
#' @param filename String - file to check
#' @param shouldExist Logical - whether you expect the file to exist; default TRUE
#' @param extension String - extension the file should have; default NA
#' @param verbose Logical - whether to print progress / results
#' @return Logical - whether file checks out
#' @details Checks whether specified file exists and returns TRUE if the file exists and was supposed to or
#' doesn't and wasn't supposed to, and FALSE otherwise.  Will also return FALSE if given an extension and the
#' file doesn't have it.
#' @author Emma Myers
#' @export



file_checks = function(filename, shouldExist=TRUE, extension=NA, verbose=FALSE) {

    # Given a filename and skip condition (either it already exists or it doesn't exist), say you're skipping if it meets the skip condition.
    # Optionally also give an extension and say you're skipping if it doesn't have that extension.

    # We'll return the value of checksOut
    checksOut = TRUE

    # Does the file exist?
    doesExist = file.exists(filename)

    # If it should and doesn't
    if ( shouldExist & !doesExist ) {
        if (verbose) { writeLines(paste("File does not exist:", filename)) }
        checksOut = FALSE
    }

    # If it does and shouldn't
    if ( doesExist & !shouldExist ) {
        if (verbose) { writeLines(paste("File already exists:", filename)) }
        checksOut = FALSE
    }

    # Are we checking for a certain extension?
    # Only check if we made it this far; otherwise the multiple messages can be confusing.
    if ( checksOut ) {
        if ( ! is.na(extension) ) {
            # If so, make it start with a dot, for consistency
            if (! substr(extension, 1, 1) == ".") {
                extension = paste('.', extension, sep='')
            }
            # Check that the last bit of the file name is .<extension>
            extLen = nchar(extension)
            if ( ! substr(filename, (nchar(filename)-extLen)+1, nchar(filename)) == extension) {
                if (verbose) { writeLines(paste("File doesn\'t have extension \'", extension, "\': ", filename, sep="")) }
                checksOut = FALSE
            }

        }
    }

    return(checksOut)

}

