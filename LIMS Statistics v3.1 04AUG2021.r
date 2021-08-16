#--------------------------------------------------------------------------------#
#   LIMS Study Statistics (Method)
#--------------------------------------------------------------------------------#

Statistics <- function(from="",to="",dataExplorer=FALSE){
    openFile(dataExplorer)
    if (dataExplorer){
        convExplorer()
    } else {
        convLIMS()
    }
    filterByDate(from, to)
    studyStatistics()
    bioStatistics()
}

#--------------------------------------------------------------------------------#
#   Load Libraries
#--------------------------------------------------------------------------------#

library(stringr)
library(dplyr)

#--------------------------------------------------------------------------------#
#   Open LIMS File
#--------------------------------------------------------------------------------#

openFile <- function(dataExplorer){
    lims <<- read.table(
        gsub("\\\\", "/", 
            choose.files(
                default = "M:/04 Studies/4.1 Active Studies/)",
                caption = "Select export_to_sponsor", 
                multi = F
            )
        ),  
        header = T,
        sep = ifelse(dataExplorer, ";", "|")
    )
}

#--------------------------------------------------------------------------------
#   Convert LIMS to Table
#--------------------------------------------------------------------------------

convLIMS <- function(){
    resList <- c("POSITIVE", "NEGATIVE", "P", "N", "ND")

    # Collect All Comments
    com <<- subset(lims, 
        grepl("[A-z]", FORMATTED_ENTRY) & 
        !FORMATTED_ENTRY %in% resList
    )
    com <<- rename(com, "COMMENT" = "FORMATTED_ENTRY")
    com <<- com[, c("REQUISITION", "ANALYSIS", "COMMENT")]

    # Collect All Results
    res <<- subset(lims, 
        grepl("[.|0-9|]", FORMATTED_ENTRY) &
        !grepl("[-/]", FORMATTED_ENTRY) |
        FORMATTED_ENTRY %in% resList
    )
    res <<- rename(res, "RESULT" = "FORMATTED_ENTRY")
    res <<- res[, c("PROJECT", "REQUISITION", "ANALYSIS", "RESULT")]

    # Collect All Dates
    dates <<- subset(lims, 
        grepl("[-|/]", FORMATTED_ENTRY) & 
        !grepl("is below", FORMATTED_ENTRY)
    )
    dates <<- rename(dates, "ANADATE" = "FORMATTED_ENTRY")
    dates <<- dates[,c("REQUISITION", "ANALYSIS", "ANADATE")]

    # Rejoin LIMS Data
    lims <<- left_join(res, dates, by=c("REQUISITION", "ANALYSIS"))
    lims <<- left_join(lims, com, by=c("REQUISITION", "ANALYSIS"))
}

#--------------------------------------------------------------------------------
#   Convert Data Explorer to Table
#--------------------------------------------------------------------------------

convExplorer <- function(){
    # Collect Comments
    com <<- subset(lims, NAME=="Comment")
    com <<- rename(com, c(
        "PROJECT" = "X_PROJECT",
        "COMMENT" = "ENTRY", 
        "REQUISITION" = "TEXT_ID"
    ))
    com <<- com[,c("ANALYSIS", "COMMENT", "REQUISITION")]

    # Collect Results
    res <<- subset(lims, NAME=="concentration")
    res <<- rename(res, c(
        "PROJECT" = "X_PROJECT", 
        "RESULT" = "ENTRY", 
        "REQUISITION" = "TEXT_ID"
    ))
    res <<- res[,c("PROJECT", "ANALYSIS", "RESULT", "REQUISITION")]

    # Collect Analysis Date
    dates <<- subset(lims, NAME=="analysis date")
    dates <<- rename(dates, c(
        "PROJECT" = "X_PROJECT",
        "ANADATE" = "ENTRY", 
        "REQUISITION" = "TEXT_ID"
    ))
    dates <<- dates[,c("ANALYSIS", "ANADATE", "REQUISITION")]
    
    # Rejoin LIMS Data
    lims <<- left_join(res, com, by=c("REQUISITION", "ANALYSIS"))
    lims <<- left_join(lims, dates, by=c("REQUISITION", "ANALYSIS"))
}

#--------------------------------------------------------------------------------#
#   Filter LIMS by Date
#--------------------------------------------------------------------------------#

filterByDate <- function(from, to){
    # Filter Data by Dates
    lims$DATE <<- format(as.Date(lims$ANADATE, format="%d-%m-%Y"), format="%Y-%m-%d")
    
    if (from!=""){
        from <- format(as.Date(from,format="%d%b%y"), format="%Y-%m-%d")
        lims <<- subset(lims, DATE>=as.character(from))
    }

    if (to!=""){
        to <- format(as.Date(to,format="%d%b%y"), format="%Y-%m-%d")
        lims <<- subset(lims, DATE<=as.character(to))
    }
}

#--------------------------------------------------------------------------------#
#   Study Statistics
#--------------------------------------------------------------------------------#

studyStatistics <- function(){
    StudyStats <<- data.frame(
        Study = lims$PROJECT[1],
        Bios = length(unique(lims$ANALYSIS)),
        N = nrow(lims),
        REQs = length(unique(lims$REQUISITION)),
        COMs = nrow(subset(lims, COMMENT!="")),
        From = format(min(as.Date(strptime(lims$ANADATE, format="%d-%m-%Y")),na.rm=T), format="%d%b%y"),
        To = format(max(as.Date(strptime(lims$ANADATE, format="%d-%m-%Y")),na.rm=T), format="%d%b%y")
    )
    cat("\n------------------------------------------------------------\n")
    print(StudyStats)
    cat("------------------------------------------------------------\n")
}

#--------------------------------------------------------------------------------#
#   Biomarker Statistics
#--------------------------------------------------------------------------------#

bioStatistics <- function(dataExplorer){
    BioStats <<- setNames(
        data.frame(matrix(nrow=0,ncol=8)),
        c("Bio", "N", "REQs", "COMs", "LLOQ", "Hemo", "From", "To")
    )
    bios <- unique(lims$ANALYSIS)

    for (bio in bios){
        subBio <- subset(lims, ANALYSIS==bio)
        BioStats[nrow(BioStats)+1,] <<- c(
            bio,
            nrow(subBio),
            length(unique(subBio$REQUISITION)),
            nrow(subset(subBio, COMMENT!="")),
            nrow(subset(subBio, tolower(COMMENT)=="below lower limit of quantification" & tolower(COMMENT)=="below_lloq")),
            nrow(subset(subBio, tolower(COMMENT)=="hemolysis" | tolower(COMMENT)=="severe hemolysis")),
            format(min(as.Date(strptime(subBio$ANADATE, format="%d-%m-%Y")),na.rm=T), format="%d%b%y"),
            format(max(as.Date(strptime(subBio$ANADATE, format="%d-%m-%Y")),na.rm=T), format="%d%b%y")
        )
    }
    print(BioStats)
    cat("------------------------------------------------------------\n")
}

#--------------------------------------------------------------------------------#
#   Run Method
#--------------------------------------------------------------------------------#

Statistics()
