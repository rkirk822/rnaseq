
# rename_files.R
# 
# This is for a very specific situation where you want to re-name a bunch of files and
# you have a csv file with two columns with the old and new filenames.  You can tell it
# to assume there's a suffix to the filenames, but it has to be the same for all files
# (both old and new names).

rename_files = function(nameFile, fileSuffix='', header=FALSE, comment.char='#') {

   # Get the old and new file prefixes and append suffix
    allNames = read.csv(nameFile, comment.char=comment.char, header=header)
    paster=function(x){paste(x, fileSuffix, sep='')}
    allNames = sapply(allNames, paster)
    
    for (i in 1:dim(allNames)[1]) {
        
        oldName = allNames[i,1]
        newName = allNames[i,2]
        
        # Check that old file exists
        if ( !file.exists(oldName) ) {
            writeLines(paste('File', oldName, 'does not exist; skipping.'))
            next
        }
    
        # Don't overwrite anything
        if ( file.exists(newName) ) {
            writeLines(paste('File', oldName, 'will not be re-named, as file', newName, 'already exists.'))
            next
        }
    
        writeLines(paste('Re-naming', oldName, 'to', newName))
        file.rename(oldName, newName)
    
    }

}

