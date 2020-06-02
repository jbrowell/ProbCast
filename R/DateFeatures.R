#' Load data of UK Holidays
#'
#' Function to read official data of UK public holidays.
#'
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' @param folder Folder location of .ics file(s) of UK holidays
#' from https://www.gov.uk/bank-holidays
#' @return Table of UK Holidays
#' @keywords Date Features
#' @import data.table
#' @export
read_holiday_ics <- function(folder=getwd()){

  for(file in list.files(path = folder,pattern = ".ics",full.names = T)){

    x <- readLines(file,encoding =  "UTF-8")
    x <- gsub(x,pattern = "Andrew.+s",replacement = "Andrew's")
    x <- gsub(x,pattern = "Year.+s",replacement = "Year's")
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
#' Add calendar-based columns to a \code{data.table} with a date/time
#' column, e.g. "hour-of-day" or holidays (based on \code{read_holiday_ics()}).
#'
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' @param dt \code{data.table} to have date/time feature columns added
#' @param datetimecol column name of \code{POSIX} time stamps
#' @param UKHolidays optional table of holidays from \code{read_holiday_ics}
#' @return Table of UK Holidays
#' @keywords Date Features
#' @import data.table
#' @export
add_calendar_variables <- function(dt,datetimecol,UKHolidays=NULL){
  
  if(dt[,attributes(get(datetimecol))$tz]=="UTC"){warning("Time zone is UTC. Local time may be prefereable for clock_hour.")}
  
  # time-since beginning of dataset
  dt[,t:=as.numeric(get(datetimecol)-min(get(datetimecol)))]
  dt[,t:=t/max(t)]
  
  
  # Day-of-week
  dt[,dow:= as.factor(format(get(datetimecol),"%a"))]
  
  # Day-of-year
  dt[,doy:= yday(get(datetimecol))]
  
  # Clock Hour
  dt[,clock_hour:=hour(get(datetimecol))+minute(get(datetimecol))/60]
  
  # UK Holidays
  if(!is.null(UKHolidays)){
    
    # UKHolidays <- read_holiday_ics()
    UKHolidays[,Date:=as.Date(Date)]
    UKHolidays <- UKHolidays[Date>=dt[,min(get(datetimecol))] & Date<=dt[,max(get(datetimecol))],]
    
    
    dt[,Date:=as.Date(get(datetimecol))]
    dt[,c("type","hol_EW","hol_Sc"):=list("N","N","N")]
    for(i in 1:nrow(UKHolidays)){
      dt[Date==UKHolidays[i,Date],c("hol_EW","hol_Sc"):=list(UKHolidays[i,`england-and-wales`],
                                                            UKHolidays[i,scotland])]
    }
    dt[is.na(hol_EW),hol_EW:="N"]; dt[is.na(hol_Sc),hol_Sc:="N"]
    dt[hol_EW!="N" & hol_Sc!="N",type:="UK"]
    dt[hol_EW=="N" & hol_Sc!="N",type:="Sc"]
    dt[hol_EW!="N" & hol_Sc=="N",type:="EW"]
    
    ## Correct for substitute days
    
    # Christmas Day
    dt[hol_EW=="Christmas Day" & (day(Date)!=25 | month(Date)!=12),hol_EW:="Substitute"]
    dt[hol_Sc=="Christmas Day" & (day(Date)!=25 | month(Date)!=12),hol_Sc:="Substitute"]
    dt[day(Date)==25 & month(Date)==12,c("hol_EW","hol_Sc"):="Christmas Day"]
    
    # Boxing Day
    dt[hol_EW=="Boxing Day" & (day(Date)!=26 | month(Date)!=12),hol_EW:="Substitute"]
    dt[hol_Sc=="Boxing Day" & (day(Date)!=26 | month(Date)!=12),hol_Sc:="Substitute"]
    dt[day(Date)==26 & month(Date)==12,c("hol_EW","hol_Sc"):="Boxing Day"]
    
    # New Year's Day
    dt[hol_EW=="New Year's Day" & (day(Date)!=1 | month(Date)!=1),hol_EW:="Substitute"]
    dt[hol_Sc=="New Year's Day" & (day(Date)!=1 | month(Date)!=1),hol_Sc:="Substitute"]
    dt[day(Date)==1 & month(Date)==1,c("hol_EW","hol_Sc"):="New Year's Day"]
    
    # 2nd January [Scotland Only]
    dt[hol_Sc=="2nd January" & (day(Date)!=2 | month(Date)!=1),hol_Sc:="Substitute"]
    dt[day(Date)==2 & month(Date)==1,hol_Sc:="2nd January"]
    
    # St Andrew's Day [Scotland Only]
    # dt[hol_Sc=="St Andrew's Day" & (day(Date)!=30 | month(Date)!=11),hol_Sc:="Substitute"]
    # dt[day(Date)==30 & month(Date)==11,hol_Sc:="2nd January"]
    
    ## Not enough examples... NYD is similar to substitute days
    dt[hol_Sc=="Substitute",hol_Sc:="New Year's Day"]
    dt[hol_EW=="Substitute",hol_EW:="New Year's Day"]
  
    ## Add bridging days? <<<<
    

  }
  
  return(dt)
}

