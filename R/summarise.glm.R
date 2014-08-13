summarise.glm <- function(
  lstModels,
  outfunc=exp,
  writetab=TRUE,
  file="modsum.csv",
  sigdigits=3,
  transpose=FALSE) {
  
  # Figure out the number of models
  nomods <- length(lstModels) 
  
  # Make a vector of all coefficients
  allCoeffs <- c()
  for (i in 1:nomods) {
    
    # Select current model results
    mod <- lstModels[[i]]
    
    # Get a list of variables
    vars <- names(mod$coefficients)
    novars <- length(vars)
    
    # Go through each variabel and add it if its not already in the list
    for (j in 1:novars) {
      
      # Get the variable name
      curname <- vars[j]
      
      # Test for the presence of the variable in the master list
      var_present <- (curname %in% allCoeffs)
      
      # If not in the list add it
      if (!(var_present)) allCoeffs <- c(allCoeffs,curname)
      
      # Close the for loop for j
    }
    
    # Close the for loop for i	
  }
  
  # Define the data structures used to extract the information from the models 	
  noCoeffs <- length(allCoeffs)
  matPointEst <- matrix(NA,nrow=nomods,ncol=noCoeffs,dimnames=list(1:nomods,allCoeffs))
  matLB <- matrix(NA,nrow=nomods,ncol=noCoeffs,dimnames=list(1:nomods,allCoeffs))
  matUB <- matrix(NA,nrow=nomods,ncol=noCoeffs,dimnames=list(1:nomods,allCoeffs))
  vecAIC <- vector(mode="numeric",length=nomods)
  vecDEX <- vector(mode="numeric",length=nomods)
  
  # Loop back though the models and the coeffciients to populate the data structures
  for (i in 1:nomods) {
    
    # Select current model results
    mod <- lstModels[[i]]
    cis <- confint.default(mod)
    
    # Get a list of variables
    vars <- names(mod$coefficients)
    novars <- length(vars)
    
    # Record the AIC
    vecAIC[i] <- mod$aic
    vecDEX[i] <- (1-mod$deviance/mod$null.deviance)*100
    
    # Go through each variabel and add it if its not already in the list
    for (j in 1:novars) {
      
      # Get the variable name
      curname <- vars[j]
      
      # Extract the point estimate and confidence intervals for the parameters
      matPointEst[i,curname] 	<- mod$coefficients[curname]  
      matLB[i,curname] 		<- cis[curname,1]  
      matUB[i,curname] 		<- cis[curname,2]  
      
      # Close the for loop for j
    }
    
    # Close the for loop for i	
  }
  
  # If selected, write a nicely formatted csv table for the parameters and models
  if (writetab) {
    
    if (transpose) {
      
      # Declare the output string
      strTable <- ""
      
      # Put in the first header row
      strTable <- paste(strTable,"Parameter",sep="")
      for (i in 1:noCoeffs) strTable <- paste(strTable,",",allCoeffs[i],",",allCoeffs[i],sep="")
      strTable <- paste(strTable,",AIC,DEX\n",sep="")
      
      # Put in the second header row
      strTable <- paste(strTable,"Model",sep="")
      for (i in 1:noCoeffs) strTable <- paste(strTable,",PE,CI",sep="")
      strTable <- paste(strTable,",AIC,DEX\n",sep="")
      
      # Output individual model lines, starting with coefficient loop
      for (i in 1:nomods) {
        
        # Pull the name of the current coefficient
        # curname <- allCoeffs[i]
        
        # Put in the name of the coefficient
        strTable <- paste(strTable,i,sep="")
        
        # Cycle through the tables looking at the different models
        for (j in 1:noCoeffs) {
          
          # Itentify the current coefficient
          curname <- allCoeffs[j]
          
          # Put in the point estimates and confidence intervals for each parameter / model combination
          curPE <- signif(outfunc(matPointEst[i,curname]),digits=sigdigits)
          curLB <- signif(outfunc(matLB[i,curname]),digits=sigdigits)
          curUB <- signif(outfunc(matUB[i,curname]),digits=sigdigits)
          
          # Paste in the parameter values and the confidence intervals
          if (is.na(curPE)) {
            
            # Put in the entry for NA results
            strTable <- paste(strTable,",","-",",","-",sep="")
            
          } else {
            
            # Put in the entry for non NA results
            strTable <- paste(strTable,",",curPE,",","(",curLB,"--",curUB,")",sep="")
            
          }
          
          # End j loop for coefficients
        }
        
        # Add the AIC at the end of the line, with a return
        mod <- lstModels[[i]]
        curAIC <- round(mod$aic,digits=1)
        curDEX <- round((1-mod$deviance/mod$null.deviance)*100,digits=1)
        strTable <- paste(strTable,",",curAIC,",",curDEX,"\n",sep="")
        
        # End the i for loop for models
      }
      
      # End the if clause for transpose
    } else {
      
      # Declare the output string
      strTable <- ""
      
      # Put in the first header row
      strTable <- paste(strTable,",Model 1",sep="")
      if (nomods>1) for (i in 2:nomods) strTable <- paste(strTable,",,Model ",i,sep="")
      strTable <- paste(strTable,"\n",sep="")
      
      # Put in the second header row
      if (nomods>1) for (i in 1:nomods) {
          strTable <- paste(strTable,",Estimate,(95% CI)",sep="")
        }
      strTable <- paste(strTable,"\n",sep="")
      
      # Output individual coefficient lines, starting with coefficient loop
      for (i in 1:noCoeffs) {
        
        # Pull the name of the current coefficient
        curname <- allCoeffs[i]
        
        # Put in the name of the coefficient
        strTable <- paste(strTable,curname,sep="")
        
        # Cycle through the tables looking at the different models
        for (j in 1:nomods) {
          
          # Put in the point estimates and confidence intervals for each 
          # parameter / model combination
          curPE <- signif(outfunc(matPointEst[j,curname]),digits=sigdigits)
          curLB <- signif(outfunc(matLB[j,curname]),digits=sigdigits)
          curUB <- signif(outfunc(matUB[j,curname]),digits=sigdigits)
          
          # Paste in the parameter values and the confidence intervals
          if (is.na(curPE)) {
            
            # Put in the entry for NA results
            strTable <- paste(strTable,",","-",",","-",sep="")
            
          } else {
            
            # Put in the entry for non NA results
            strTable <- paste(strTable,",",curPE,",","(",curLB,"--",curUB,")",sep="")
            
          }
          
          # End model for loop
        }
        
        # Return at the end of the line
        strTable <- paste(strTable,"\n",sep="")
        
        # End for for coeffs	
      }
      
      
      # Write the row name for the AICs	
      strTable <- paste(strTable,"AIC",sep="")
      
      # Start the for loop for the AICs
      for (i in 1:nomods) {
        
        # Get the current model
        mod <- lstModels[[i]]
        
        # Format the AIC for the current model
        curAIC <- round(mod$aic,digits=1)
        
        # Write the value and the space
        strTable <- paste(strTable,",",curAIC,",",sep="")
      }
      
      # Return at the end of the AIC line
      strTable <- paste(strTable,"\n",sep="")
      
      # Write the row name for the DEXs	
      strTable <- paste(strTable,"DEX",sep="")
      
      # Start the for loop for the DEX
      for (i in 1:nomods) {
        
        # Get the current model
        mod <- lstModels[[i]]
        
        # Format the AIC for the current model
        curDEX <- round((1-mod$deviance/mod$null.deviance)*100,digits=1)
        
        # Write the value and the space
        strTable <- paste(strTable,",",curDEX,",",sep="")
      }
      
      # Return at the end of the DEX line
      strTable <- paste(strTable,"\n",sep="")
      
      # End else statement for transpose of 
    }
    
    # Write the string to the selected file
    cat(strTable,file=file)
    
  }
  
  # Return 
  data.frame(pe=matPointEst,lb=matLB,ub=matUB,aic=vecAIC,dex=vecDEX)
  
}
