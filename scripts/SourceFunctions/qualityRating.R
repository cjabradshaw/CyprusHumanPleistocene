qualityRating <- function(x) {
  
  preQuality <- character(nrow(x))
  Reason <- character(nrow(x))
  
  for (i in 1:nrow(x)) {
    a <- x[i, ]
    
    # If there is no lab code
    if (a$AgeID == "?"| 
        a$AgeID =="na" | 
        a$AgeID =="not reported" | 
        a$AgeID =="Not given") {
      preQuality[i] <- "C"
      Reason[i] <- "Missing lab code"
    }
    
    # If the dating technique is not known
    else if (a$DatingTechnique == "na") {
      preQuality[i] <- "C"
      Reason[i] <- "Dating technique unknown"
    }
    
    # If it is an age range or presence only data
    else if (a$AgeType == "range") {
      preQuality[i] <- "C"
      Reason[i] <- "Date range"
    }
    
    else if (a$AgeType == "presence") {
      preQuality[i] <- "C"
      Reason[i] <- "Presence only"
    }
    
    
    # If it is a radiocarbon age
    else if (grepl("Radiocarbon", a$DatingTechnique)) {
      
      
      
      # If dating was made on a mixture of bones
      if (grepl("mixed bones", a$DatedMaterial, ignore.case = TRUE) |
          grepl("mixed bones", a$DatedRemain, ignore.case = TRUE) |
          grepl("bulk bones", a$DatedMaterial, ignore.case = TRUE) |
          grepl("bulk bones", a$DatedRemain, ignore.case = TRUE) ) {
        preQuality[i] <- "C"
        Reason[i] <- "Mixed/bulk bones were dated"
      }
      
      # If bulk soil organics were dated
      else if (grepl("bulk soil organic", a$DatedMaterial, ignore.case = TRUE) |
               grepl("bulk soil organic", a$DatedRemain, ignore.case = TRUE)) {
        preQuality[i] <- "C"
        Reason[i] <- "Bulk soil organics were dated"
      }
      
      # If inorganic calcite was dated
      else if (grepl("Speleothem", a$DatedMaterial, ignore.case = TRUE) |
               grepl("Speleothem", a$DatedRemain, ignore.case = TRUE) |
               grepl("Flowstone", a$DatedRemain, ignore.case = TRUE) |
               a$DatedRemain == "Sediment" & 
               a$DatedMaterial == "Carbonate") {
        preQuality[i] <- "C"
        Reason[i] <- "Inorganic calcite dated"
      }
      
      # If bone or teeth were dated
      else if (grepl("Bone", a$DatedRemain, ignore.case = TRUE) |
               grepl("Tooth", a$DatedRemain, ignore.case = TRUE) |
               grepl("Teeth", a$DatedRemain, ignore.case = TRUE)) {
        
        # If single aminoacid was dated
        if (a$DatedMaterial == "Single aminoacid") {
          
          # If extraction problems were reported
          if (a$C14_ExtractionProblem == "Yes") {
            preQuality[i] <- "B"
            Reason[i] <- "Reported extraction problems"
          }
          
          #if pretreatment was adequate
          else if (grepl("XAD-2", a$C14_Pretreatment, ignore.case = TRUE) |
                   grepl("Ninhydrin", a$C14_Pretreatment, ignore.case = TRUE) | 
                   grepl("ultrafiltration", a$C14_Pretreatment, ignore.case = TRUE) |
                   grepl("C-AF", a$C14_Pretreatment, ignore.case = TRUE)) {
            preQuality[i] <- "m*"
            Reason[i] <- "Individual aminoacids dated with adequate pretreatment"
          }
          
          else {
            preQuality[i] <- "B"
            Reason[i] <- "Inadequate pretreatment"
          }
        }
        
        # If collagen was dated
        else if (a$DatedMaterial == "Collagen") {
          
          # If extraction problems were reported
          if (a$C14_ExtractionProblem == "Yes") {
            preQuality[i] <- "B"
            Reason[i] <- "Reported problems with collagen extraction"
          }
          
          # If collagen preservation is unknown
          else if (a$C14_CNRatioValue == "Not reported" | 
                   a$C14_CNRatioValue == "na") {
            preQuality[i] <- "B"
            Reason[i] <- "Unknown collagen preservation" 
          }
          
          #if pretreatment was adequate
          else if (grepl("XAD-2", a$C14_Pretreatment, ignore.case = TRUE) |
                   grepl("Ninhydrin", a$C14_Pretreatment, ignore.case = TRUE) | 
                   grepl("ultrafiltration", a$C14_Pretreatment, ignore.case = TRUE) |
                   grepl("C-AF", a$C14_Pretreatment, ignore.case = TRUE)){
            
            #if CN ratio was reported
            if (a$C14_CNRatioValue>2 && a$C14_CNRatioValue <4) {
              
              #if %N was reported
              if (a$C14_NPercentage > 0){
                preQuality[i] <- "m*"
                Reason[i] <- "All requirements met"
              }
              
              else {
                preQuality[i] <- "B"
                Reason[i] <- "%N not reported"
              }
            }
            
            else {
              preQuality[i] <- "B"
              Reason[i] <- "CN ratio not reported"
            }
          }
          
          #if pretreatment is inadequate
          else {
            preQuality[i] <- "B"
            Reason[i] <- paste("Pretreatment inadequate (", a$C14_Pretreatment, ")")
          }
          
        }#closes collagen section (radiocarbon)
        
        # If tooth enamel was dated
        else if (a$DatedMaterial == "Enamel") {
          preQuality[i] <- "C"
          Reason[i] <- "Tooth enamel was dated"
        }
        
        #If apatite was dated
        else if (a$DatedMaterial == "Apatite") {
          preQuality[i] <- "C"
          Reason[i] <- "Bone apatite was dated"
        }
        
        # If bone was given as dated material but no further specification
        else {
          preQuality[i] <- "C"
          Reason[i] <- "Inadequate material was dated"
        }
        
      }#closes bone and teeth section (radiocarbon)
      
      #If wood or seeds were dated
      else if (a$DatedRemain == "wood" | 
               a$DatedRemain == "seed") {
        
        if (a$DatedMaterial == "alpha-cellulose") {
          preQuality[i] <- "m"
          Reason[i] <- "All requirements met"
        }
        
        else {
          preQuality[i] <- "C"
          Reason[i] <- "Unknown molecule used for dating"
        }
        
      }#closes woods and seeds section (radiocarbon)
      
      #If eggshell was dated
      else if (grepl("eggshell", a$DatedRemain, ignore.case = TRUE)) {
        
        if (a$DatedMaterial == "Organic fraction") {
          preQuality[i] <- "C"
          Reason[i] <- "Organic fraction was dated"
        }
        
        else if (a$DatedMaterial == "Carbonate") {
          
          if (a$C14_Pretreatment == "Grinding and acid etching") {
            preQuality[i] <- "m*"
            Reason[i] <- "All requirements met"
          }
          
          else {
            preQuality[i] <- "B"
            Reason[i] <- paste("Inadequate pretreatment to remove seconday carbonate (", a$C14_Pretreatment, ")")
          }
          
        }
        
        else {
          preQuality[i] <- "C"
          Reason[i] <- "Inappropriate material was dated using Radiocarbon dating"
          
        }  
        
        
      }#closes eggshell section (radiocarbon)
      
      #if shells were dated
      else if (grepl("shell", a$DatedRemain, ignore.case = TRUE) |
               grepl("coral", a$DatedRemain, ignore.case = TRUE) |
               grepl("otolith", a$DatedRemain, ignore.case = TRUE) |
               grepl("shell", a$DatedMaterial, ignore.case = TRUE) |
               grepl("coral", a$DatedMaterial, ignore.case = TRUE) |
               grepl("otolith", a$DatedMaterial, ignore.case = TRUE)){
        
        if (grepl("carbonate", a$DatedMaterial, ignore.case = TRUE)) {
          
          if (a$C14_Pretreatment == "Physical cleaning and acid etching") {
            
            if(a$C14_corals/shells_X-ray_diffraction_recrystallisation == "yes") {
              preQuality[i] <- "m"
              Reason[i] <- "All requirements met"
            }
            
            else {
              preQuality[i] <- "B"
              Reason[i] <- "Shell, coral or otolith was dated without X-ray_diffraction"
            }
            
          }
          
          else {
            preQuality[i] <- "B"
            Reason[i] <- paste("Shell, coral or otolith carbonate was dated with insufficient pretreatment:", a$C14_Pretreatment)
          }
        }
        
        else {
          preQuality[i] <- "B"
          Reason[i] <- "Shell, coral or otolith was dated without further specification"
        }
        
      }#closes shell section (radiocarbon)
      
      #if gut contents of coprolites were dated
      else if (grepl("gut contents", a$DatedRemain, ignore.case = TRUE) |
               grepl("coprolite", a$DatedRemain, ignore.case = TRUE)) {
        
        if (a$DatedMaterial == "alpha-cellulose") {
          preQuality[i] <- "m"
          Reason[i] <- "All requirements met"
        }
        
        else {
          preQuality[i] <- "C"
          Reason[i] <- "Unknown molecule used for dating"
        }
      }#closes gut content and coprolite section
      
      
      #If charcoal was dated
      else if (grepl("Charcoal", a$DatedMaterial, ignore.case = TRUE)) {
        
        #If no info about pretreatment is given
        if (a$C14_Pretreatment == "Not reported" | 
            a$C14_Pretreatment == "not reported" | 
            a$C14_Pretreatment == "na"){
          preQuality[i] <- "C"
          Reason[i] <- "Missing info about pretreatment"
        }
        
        else if (grepl("ABOX", a$C14_Pretreatment, ignore.case = TRUE) |
                 grepl("chlorate oxidation", a$C14_Pretreatment, ignore.case = TRUE)) {
          preQuality[i] <- "m"
          Reason[i] <- "Adequate pretreatment"
        }        
        
        else {
          preQuality[i] <- "B"
          Reason[i] <- paste("Pretreatment inadequate (", a$C14_Pretreatment, ")")
        }
        
      }#closes charcoal section (radiocarbon)
      
      # If no info about dated molecule was provided
      else if (a$DatedMaterial == "Not reported" | 
               a$DatedMaterial == "") {
        preQuality[i] <- "C"
        Reason[i] <- "No info about dated molecule"
      }
      
      # If anything else was dated
      else {
        preQuality[i] <- "C"
        Reason[i] <- "Inappropriate material was dated using Radiocarbon dating"
      }
      
      
    } #closes radiocarbon section
    
    
    #if it is an AAR age
    else if (grepl("AAR", a$DatingTechnique)) {
      
      #if bone or tooth was dated
      if (#grepl("bone", a$DatedMaterial, ignore.case = TRUE) |
        grepl("bone", a$DatedRemain, ignore.case = TRUE) |
        #grepl("tooth", a$DatedMaterial, ignore.case = TRUE) |
        grepl("tooth", a$DatedRemain, ignore.case = TRUE) |
        #grepl("teeth", a$DatedMaterial, ignore.case = TRUE) |
        grepl("teeth", a$DatedRemain, ignore.case = TRUE)) {
        preQuality[i] <- "C"
        Reason[i] <- "Bone or teeth were dated using AAR"
      }
      
      else if (grepl("eggshell", a$DatedRemain, ignore.case = TRUE) |
               grepl("otolith", a$DatedRemain, ignore.case = TRUE)){
        
        if(grepl("burnt", a$AAR_thermal_history, ignore.case = TRUE) |
           grepl("unknown", a$AAR_thermal_history, ignore.case = TRUE)) {
          preQuality[i] <- "C"
          Reason[i] <- paste("thermal history (", a$AAR_thermal_history, ")")
        } 
        
        else if(grepl("yes", a$AAR_demonstrated_CS_behaviour, ignore.case = TRUE)){
          
          if(grepl("yes", a$AAR_replication_with_low_uncertainty, ignore.case = TRUE)){
            
            if(a$Association == "direct" | a$Association == "Direct"){
              
              if(grepl("yes", a$AAR_reliable_calibration, ignore.case = TRUE)){
                preQuality[i] <- "m*"
                Reason[i] <- "All requirements met"
              }
              
              else{
                preQuality[i] <- "m"
                Reason[i] <- "No independent calibration"
              }
            }
            
            else{
              preQuality[i] <- "C"
              Reason[i] <- "AAR indirect age"
            }
          }
          
          else{
            preQuality[i] <- "C"
            Reason[i] <- "AAR unreplicated age"
          }
        }
        
        else{
          preQuality[i] <- "C"
          Reason[i] <- "AAR on potential open system"
        }
      }#closes AAR eggshell/otholith  section
      
    }#closes AAR section
    
    
    #if it is an uranium-series age
    else if (grepl("U-TH", a$DatingTechnique, ignore.case = TRUE) |
             grepl("Th/U", a$DatingTechnique, ignore.case = TRUE)) {
      
      #if teeth were dated
      if (grepl("teeth", a$DatedMaterial, ignore.case = TRUE) |
          grepl("tooth", a$DatedMaterial, ignore.case = TRUE) |
          grepl("teeth", a$DatedRemain, ignore.case = TRUE) |
          grepl("tooth", a$DatedRemain, ignore.case = TRUE)){
        
        if(grepl("ICP-MS", a$DatingTechnique, ignore.case = TRUE) |
           grepl("TIMS", a$DatingTechnique, ignore.case = TRUE)) {
          
          if(grepl("yes", a$Uth_ClosedSystem, ignore.case = TRUE)) {
            preQuality[i] <- "m"
            Reason[i] <- "All requirements met"
          }
          
          else{
            preQuality[i] <- "B"
            Reason[i] <- "Tooth was dated without demonstrated closed-system behaviour"
          } 
        }
        
        else{
          preQuality[i] <- "C"
          Reason[i] <- "Tooth was dated without appropriate method"
        }
        
      }#closes tooth section (U-series)
      
      #if bone was dated
      else if (grepl("bone", a$DatedMaterial, ignore.case = TRUE)){
        
        if(grepl("ICP-MS", a$DatingTechnique, ignore.case = TRUE) |
           grepl("TIMS", a$DatingTechnique, ignore.case = TRUE)) {
          
          if(grepl("yes", a$Uth_ClosedSystem, ignore.case = TRUE)) {
            preQuality[i] <- "m"
            Reason[i] <- "All requirements met"
          }
          
          else{
            preQuality[i] <- "B"
            Reason[i] <- "Bone was dated without modelling"
          } 
        }
        
        else{
          preQuality[i] <- "C"
          Reason[i] <- "Bone was dated without further specification"
        }
        
      }#closes bone section (U-series)
      
      #if eggshell was dated
      else if (grepl("eggshell", a$DatedMaterial, ignore.case = TRUE)) {
        
        #if ICP-MS or TIMS was used
        if (grepl("TIMS", a$DatingTechnique, ignore.case = TRUE) |
            grepl("ICP-MS", a$DatingTechnique, ignore.case = TRUE)) {
          
          #if the pretreatment was adequate  
          if (a$Uth_Pretreatment == "Grinding and acid etching") {
            preQuality[i] <- "m*"
            Reason[i] <- "All requirements met"
          }
          
          # If no info about pretreatment is given
          else if (a$Uth_Pretreatment == "Not reported" | 
                   a$Uth_Pretreatment == "not reported" | 
                   a$Uth_Pretreatment == "na" | 
                   a$Uth_Pretreatment == "null"){
            preQuality[i] <- "C"
            Reason[i] <- "Missing info about pretreatment"
          }
          
          else {
            preQuality[i] <- "C"
            Reason[i] <- paste("Pretreatment inadequate (", a$Uth_Pretreatment, ")")
          }
        }
        
        else {
          preQuality[i] <- "C"
          Reason[i] <- "Eggshell was dated without ICP-MS or TIMS"
        }
        
      }#closes eggshell section (U-series)
      
      #if closed-system of no body remains were dated
      else if (grepl("speleothem", a$DatedMaterial, ignore.case = TRUE) |
               grepl("coral", a$DatedMaterial, ignore.case = TRUE) |
               grepl("calcite", a$DatedMaterial, ignore.case = TRUE) |
               grepl("speleothem", a$DatedRemain, ignore.case = TRUE) |
               grepl("coral", a$DatedRemain, ignore.case = TRUE) |
               grepl("calcite", a$DatedRemain, ignore.case = TRUE)){
        
        #if ICP-MS or TIMS was used
        if (grepl("TIMS", a$DatingTechnique, ignore.case = TRUE) |
            grepl("ICP-MS", a$DatingTechnique, ignore.case = TRUE)) {
          
          #if detrital correction was done
          if (grepl("yes", a$Uth_DetritalCorrection, ignore.case = TRUE )){
            preQuality[i] <- "m"
            Reason[i] <- "All requirements met" 
          }
          
          else {
            preQuality[i] <- "C"
            Reason[i] <- "Closed-system of no body remains were dated without detrital correction"
          }
        }
        
        else {
          preQuality[i] <- "C"
          Reason[i] <- "Closed-system of no body remains were dated twithout ICP-MS or TIMS"
        }
        
      }#closes closed system of no body remains section (U-series)
      
      else {
        preQuality[i] <- "C"
        Reason[i] <- "Uranium-series age on inadequate material"
      }
      
    } #closes U-series -section 
    
    
    #if it is an ESR age
    else if (grepl("ESR", a$DatingTechnique, ignore.case = TRUE)){
      
      #if enamel was dated
      if (grepl("enamel", a$DatedMaterial, ignore.case = TRUE)){
        
        #if it is combined ESR and closed-system U-series modelling
        if (grepl("CSUS-ESR", a$DatingTechnique, ignore.case = TRUE)){
          preQuality[i] <- "m*"
          Reason[i] <- "All requirements met"
        }
        
        #if EU and LU ages are model dependant with a p-value derived fro a U-series estimate
        else if (a$DatingTechnique =="US-ESR"){
          preQuality[i] <- "m"
          Reason[i] <- "All requirements met"
        }
        
        else if (a$ESR_InternalDoseRate.10.. == "yes") {
          
          if (a$ESR_GammaDoseRate == "in situ") {
            preQuality[i] <- "m*"
            Reason[i] <- "All requirements met"
          }
          
          else if (a$ESR_GammaDoseRate == "assumed") {
            preQuality[i] <- "m"
            Reason[i] <- "Gamma dose rate assumed"
          }
          
          else {
            preQuality[i] <- "C"
            Reason[i] <- "No information on gamma dose rate"
          }
        }
        
        else if (a$ESR_InternalDoseRate.10.. == "no") {
          preQuality[i] <- "B"
          Reason[i] <- "Internal dose rate > 10% of total dose rate"
        }
        
        else {
          preQuality[i] <- "B"
          Reason[i] <- "Internal dose rate not given, no U-series constraint on possible history of U-uptake"
        }
        
      }#closes enamel section (ESR)
      
      
      else {
        preQuality[i] <- "C"
        Reason[i] <- "Material other than tooth enamel was dated using ESR"
      }
      
    }#closes ESR section
    
    
    #if it is an luminescence (Luminescence dating) or thermoluminescence (TL) age
    else if (grepl("OSL", a$DatingTechnique, ignore.case = TRUE) |
             grepl("TL", a$DatingTechnique, ignore.case = TRUE)) {
      
      #if organic material was dated
      if (grepl("bone", a$DatedMaterial, ignore.case = TRUE) |
          grepl("tooth", a$DatedMaterial, ignore.case = TRUE) |
          grepl("teeth", a$DatedMaterial, ignore.case = TRUE) |
          grepl("eggshell", a$DatedMaterial, ignore.case = TRUE)) {
        preQuality[i] <- "C"
        Reason[i] <- "Organic material was dated using luminescence"
      }
      
      #if sediment was dated
      if (grepl("sediment", a$DatedRemain, ignore.case = TRUE)) {
        
        if (grepl("single grain", a$Luminescence_EquivalentDoseMeasurement, ignore.case = TRUE)) {
          
          if(grepl("adequately", a$Luminescence_BleachingStatus, ignore.case = TRUE) |
             grepl("partially", a$Luminescence_BleachingStatus, ignore.case = TRUE)) {
            
            if(grepl("yes", a$LuminescenceSingleGrain_EquivalentDoseModel, ignore.case = TRUE)) {
              preQuality[i] <- "m"
              Reason[i] <- "All requirements met"
            }  
            
            else if(grepl("adequately", a$Luminescence_BleachingStatus, ignore.case = TRUE)) {
              
              preQuality[i] <- "m"
              Reason[i] <- "All requirements met"
            }    
            
            else{
              preQuality[i] <- "C"
              Reason[i] <- "Single-grain OSL age that cannot be modelled"
            }
            
          }
          
          else {
            preQuality[i] <- "C"
            Reason[i] <- "Single-grain OSL age that is not sufficiently bleached"
          }
          
          
        }#closes Luminescence dating single-grain section  
        
        else if(grepl("single aliquot", a$Luminescence_EquivalentDoseMeasurement, ignore.case = TRUE) |
                grepl("multi aliquot", a$Luminescence_EquivalentDoseMeasurement, ignore.case = TRUE)) {
          
          if(grepl("adequately", a$Luminescence_BleachingStatus, ignore.case = TRUE)){
            
            preQuality[i] <- "m"
            Reason[i] <- "All requirements met"
          }
          
          
          else{
            preQuality[i] <- "C"
            Reason[i] <- "Age for mixed or partially bleached sediment"
          }
          
        }#closes multi-grain section (Luminescence dating)
        
        if (grepl("na", a$Luminescence_EquivalentDoseMeasurement, ignore.case = TRUE)) {
          preQuality[i] <- "C"
          Reason[i] <- "Luminescence dating method not specified"
        }
        
      }#closes sediment secion (Luminescence dating)
      
      else{
        preQuality[i] <- "C"
        Reason[i] <- "Inappropriate material dated using Luminescence dating"
      }
      
    }#closes Luminescence dating section 
    
  } ### closes for loop
  
  result <- data.frame(Species = x$Species, AgeID = x$AgeID, Association =x$Association, Association.sub.category =x$Association.sub.category, preQuality = preQuality, Reason = Reason)
  
  
}  ### closes function
