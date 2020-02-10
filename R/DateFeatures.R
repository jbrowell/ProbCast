## Read .ics of bank holidays from https://www.gov.uk/bank-holidays
#' Read UK Holidays
#'
#' @param folder Folder location of .ics file(s) of UK holidays from https://www.gov.uk/bank-holidays
#' @details Details go here...
#' @return Table of UK Holidays
#' @keywords Date Features
#' @export
read_holiday_ics <- function(folder=getwd()){

  for(file in list.files(pattern = ".ics")){

    x <- readLines(file,encoding =  "UTF-8")
    stopifnot(!any(grepl("^\\s+", x))) # disregarding value fields that have linefeeds for the sake of simplicity
    keyval <- do.call(rbind, regmatches(x, regexpr(":", x, fixed = TRUE), invert = TRUE))
    keyval <- keyval[which.max(keyval[,1]=="BEGIN" & keyval[,2]=="VEVENT"):tail(which(keyval[,1]=="END" & keyval[,2]=="VEVENT"), 1),]
    keyval <- cbind.data.frame(keyval, id=cumsum(keyval[,1]=="BEGIN" & keyval[,2]=="VEVENT"))
    df <- reshape(keyval, timevar="1", idvar="id", direction = "wide")

    temp <- data.table(Date=as.POSIXct(df$`2.DTSTART;VALUE=DATE`,format="%Y%m%d",tz="Europe/London"),
                       H=df$`2.SUMMARY`)
    setnames(temp,"H",sub('\\.ics$', '', file))
    if(exists("output")){
      output <- merge(output,temp,by="Date",all=T)
    }else{
      output <- temp
    }
  }
  return(output)
}

#' Add Calendar Features
#'
#' @param data \code{data.frame} to have date/time feature columns added
#' @param datetimecol column name of \code{POSIX} time stamps
#' @param UKHolidays optional table of holidays from \code{read_holiday_ics}
#' @details Details go here...
#' @return Table of UK Holidays
#' @keywords Date Features
#' @export
add_calendar_variables <- function(data,datetimecol,UKHolidays=NULL){

  if(attributes(data[[datetimecol]])$tz=="UTC"){warning("Time zone is UTC. Local time may be prefereable for clock_hour.")}
  
  # time-since beginning of dataset
  data$t <- as.numeric(data[[datetimecol]]-data[[datetimecol]][1])
  data$t <- data$t/max(data$t)

  # Day-of-week
  data$dow <- as.factor(format(data[[datetimecol]],"%a"))

  # Day-of-year
  data$doy <- as.numeric(format(data[[datetimecol]],"%j"))

  # Clock Hour
  data$clock_hour <- as.numeric(format(data[[datetimecol]],"%H"))+
    as.numeric(format(data[[datetimecol]],"%M"))/60

  # UK Holidays
  if(!is.null(UKHolidays)){

    # UKHolidays <- read_holiday_ics()
    UKHolidays$type <- "UK" # Generic holiday - All UK
    UKHolidays$type[is.na(UKHolidays$`england-and-wales`)] <- "Sc" # Generic holiday, Scotland
    UKHolidays$type[is.na(UKHolidays$scotland)] <- "EW" # Generic holiday, Eng & Wal
    UKHolidays$type[UKHolidays$scotland=="Christmas Day"] <- "Ch" # Christmas Day

    ## Add bridging days and holiday weekends...

    UKHolidays$`england-and-wales` <- NULL
    UKHolidays$scotland <- NULL

    data$Date <- as.Date(data[[datetimecol]])
    UKHolidays$Date <- as.Date(UKHolidays[["Date"]])

    if(!is.null(data$type)){data$type <- NULL}
    data <- merge(data,UKHolidays,by="Date",all.x = T)
    data$type[is.na(data$type)] <- "N"
  }

  return(data)
}
