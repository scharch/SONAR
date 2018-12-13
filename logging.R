
my.error.fun <- function() {
    time.stamp <- date()
    err.msg <- paste0( time.stamp, " -- DIED WITH ERROR:\n", geterrmessage())
    cat(err.msg, file="output/logs/command_history.log", append=T)
    q("no", status = 1, runLast = FALSE)
}

saveCommandLine <- function( logFile, args ) {
    file.arg.name <- "--file="
    script.name <- sub(file.arg.name, "", args[grep(file.arg.name, args)])
    script.basename <- dirname(script.name)
    args.positions <- grep("--args", args) 
        
    time.stamp  <- date()
    version     <- system(paste0("git -C ",script.basename," describe --always --dirty --tags"), intern=T)

    commandLine <- paste(c(script.name,tail(args, -1*args.positions)),collapse=' ')
    message     <- paste0("\n",time.stamp," -- SONAR ",version," run with command:\n\t",commandLine,"\n")
    cat(message,  file=logFile, append=T)
}


logSuccess <- function() {
    time.stamp <- date()
    success <- paste0( time.stamp, " -- Program finished successfully\n" )
    cat(success, file="output/logs/command_history.log", append=T)	   
}