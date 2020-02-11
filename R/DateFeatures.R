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
#' @param dt \code{data.table} to have date/time feature columns added
#' @param datetimecol column name of \code{POSIX} time stamps
#' @param UKHolidays optional table of holidays from \code{read_holiday_ics}
#' @details Details go here...
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
    
    
    # dt[hol_EW=="Christmas Day",plot(LeadTime,node,col=kfold)]
    # dt[hol_EW=="Boxing Day",plot(LeadTime,node,col=kfold)]
    # dt[hol_EW=="New Year's Day",plot(LeadTime,node,col=kfold)]
    # dt[hol_EW=="Summer bank holiday",plot(LeadTime,node,col=kfold)]
    # dt[hol_EW=="Substitute",points(LeadTime,node,col=kfold,pch=15)]
    
    ## Not enough examples... fudge fix
    dt[hol_Sc=="Substitute",hol_Sc:="New Year's Day"]
    dt[hol_EW=="Substitute",hol_EW:="New Year's Day"]
    
    ## Add bridging days? <<<<
    

  }
  
  return(dt)
}
## OLD VERSION
# add_calendar_variables <- function(data,datetimecol,UKHolidays=NULL){
# 
#   if(attributes(data[[datetimecol]])$tz=="UTC"){warning("Time zone is UTC. Local time may be prefereable for clock_hour.")}
#   
#   # time-since beginning of dataset
#   data$t <- as.numeric(data[[datetimecol]]-data[[datetimecol]][1])
#   data$t <- data$t/max(data$t)
# 
#   # Day-of-week
#   data$dow <- as.factor(format(data[[datetimecol]],"%a"))
# 
#   # Day-of-year
#   data$doy <- as.numeric(format(data[[datetimecol]],"%j"))
# 
#   # Clock Hour
#   data$clock_hour <- as.numeric(format(data[[datetimecol]],"%H"))+
#     as.numeric(format(data[[datetimecol]],"%M"))/60
# 
#   # UK Holidays
#   if(!is.null(UKHolidays)){
# 
#     # UKHolidays <- read_holiday_ics()
#     UKHolidays$type <- "UK" # Generic holiday - All UK
#     UKHolidays$type[is.na(UKHolidays$`england-and-wales`)] <- "Sc" # Generic holiday, Scotland
#     UKHolidays$type[is.na(UKHolidays$scotland)] <- "EW" # Generic holiday, Eng & Wal
#     UKHolidays$type[UKHolidays$scotland=="Christmas Day"] <- "Ch" # Christmas Day
# 
#     ## Add bridging days and holiday weekends...
# 
#     UKHolidays$`england-and-wales` <- NULL
#     UKHolidays$scotland <- NULL
# 
#     data$Date <- as.Date(data[[datetimecol]])
#     UKHolidays$Date <- as.Date(UKHolidays[["Date"]])
# 
#     if(!is.null(data$type)){data$type <- NULL}
#     data <- merge(data,UKHolidays,by="Date",all.x = T)
#     data$type[is.na(data$type)] <- "N"
#     
#     
#     
#     
#     
#     ## << THIS CODE BREAKS EVERYTHING!!! >> ##
#     ## Add actual date for fixed-date holidays where substitude days are used
#     # temp1 <- melt(UKHolidays[,.(Date,`england-and-wales`,scotland)],id.vars = 1,value.name = "holiday",variable.name = "region")
#     # 
#     # temp1[holiday=="Christmas Day",c("fixed_month","fixed_day"):=list(12,25)]
#     # temp1[holiday=="Boxing Day",c("fixed_month","fixed_day"):=list(12,26)]
#     # temp1[holiday=="New Year's Day",c("fixed_month","fixed_day"):=list(1,1)]
#     # temp1[holiday=="2nd January",c("fixed_month","fixed_day"):=list(1,2)]
#     # temp1[holiday=="St Andrew's Day",c("fixed_month","fixed_day"):=list(11,30)]
#     # 
#     # temp_hol <- temp1[!is.na(fixed_day) & (month(Date)!=fixed_month | day(Date)!=fixed_day),]
#     # temp_hol[,holiday:="Substitute"]
#     # temp_fixed <- temp1[!is.na(fixed_day) & (month(Date)!=fixed_month | day(Date)!=fixed_day),]
#     # temp_fixed <- temp_fixed[!is.na(fixed_day) & (month(Date)!=fixed_month | day(Date)!=fixed_day),
#     #                          Date:=as.Date(paste0(year(Date),"-",fixed_month,"-",fixed_day),tz="UTC")]
#     # 
#     # temp1 <- rbind(temp1[!(!is.na(fixed_day) & (month(Date)!=fixed_month | day(Date)!=fixed_day)),],
#     #                temp_fixed,temp_hol)
#     # temp1 <- dcast(temp1[,.(Date,region,holiday)],Date~region,value.var = "holiday")
#     # 
#     # setnames(temp1,c("england-and-wales","scotland"),c("hol_EW","hol_Scot"))
#     # data <- merge(data,temp1,by="Date",all.x = T)
#   }
# 
#   return(data)
# }
