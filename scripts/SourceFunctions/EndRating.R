EndRating<- function(x) {
  
  #preQuality <- character(nrow(x))
  Quality <- character(nrow(x))
  QualityReason <- character(nrow(x))
  
  for (i in 1:nrow(x)) {
    a <- x[i, ]
    
    
    if(a$preQuality =="m*"){
      if(a$Association =="Direct" | a$Association =="direct"){
        Quality[i] <- "A*"
        QualityReason[i] <- "Direct age"
      }
      else if(a$Association =="Indirect" | a$Association =="indirect"){
        if (a$Association.sub.category == "Yes" | a$Association.sub.category == "yes"){
          Quality[i] <- "A"
          QualityReason[i] <- "Indirect age"
        }
        else if(a$Association.sub.category == "Uncertain" | a$Association.sub.category == "uncertain"){
          Quality[i] <- "B"
          QualityReason[i] <- "Uncertain stratographical association"
        }
        else if(a$Association.sub.category == "No" | a$Association.sub.category == "no"){
          Quality[i] <- "C"
          QualityReason[i] <- "No stratographical association"
        }
        else {
          Quality[i] <- "NA"
          QualityReason[i] <- "Association sub category not clear"
        }
      }
      else {
        Quality[i] <- "NA"
        QualityReason[i] <- "Not clear if age is direct or indirect"
      }
    }
    else if (a$preQuality =="m"){
      if(a$Association =="Direct" | a$Association =="direct"){
        Quality[i] <- "A"
        QualityReason[i] <- "Direct age"
      }
      else if (a$Association =="Indirect" | a$Association =="indirect"){ 
        if (a$Association.sub.category == "Yes" | a$Association.sub.category == "yes"){
          Quality[i] <- "A"
          QualityReason[i] <- "Indirect age"
        }
        else if(a$Association.sub.category == "Uncertain" | a$Association.sub.category == "uncertain"){
          Quality[i] <- "B"
          QualityReason[i] <- "Uncertain stratographical association"
        }
        else if(a$Association.sub.category == "No" | a$Association.sub.category == "no"){
          Quality[i] <- "C"
          QualityReason[i] <- "No stratographical association"
        } 
        else {
          Quality[i] <- "NA" 
          QualityReason[i] <- "Association sub category not clear"
        }
      }
    }
    else if (a$preQuality =="B"){
      Quality[i] <- "B"
      #  QualityReason[i] <- paste(a$Reason)
    }
    else if (a$preQuality =="C"){
      Quality[i] <- "C" 
      QualityReason[i] <- paste(a$Reason)
    }
    else {
      Quality[i] <- "NA"
      QualityReason[i] <- paste(a$Reason)
    }
  }#closes for loop
  
  rating <- data.frame(Species = x$Species, AgeID = x$AgeID, Reason = x$Reason, preQuality = x$preQuality, Association =x$Association, Association.sub.category =x$Association.sub.category, Quality = Quality, QualityReason = QualityReason)
  
}
