renamed <- function(ff) {
  
  # Extract parameters and markers name
  names=unname(ff@parameters@data$name)
  mkr=unname(ff@parameters@data$desc)
  
  # scatter parameters
  Geom = c("FSC-A","FSC-H","SSC-A","SSC-H","FSC.A","FSC.H","SSC.A","SSC.H","Time") 
  
  # markers withouth CD in their name
  NoCD = c("TDT","HLADR","LAMBDA","KAPPA","IGM","TDTc","NG2","KORSA","smCD3","cyMPO","cyTDT","cyCD79","PD1",
           "cMPO","cCD79a","cyCD3","cyIGM","IREM2","NuTdT","TCRab","TCRgd",
           "MPO","TSLPR","MPX","IGG","cMPX","cKAPPA","cLAMBDA","HLA","IG","X")
  
  # markers to delete
  TOdel = c("CD-","CD","CD--","CD---","CD----","X","x")
  
  for (i in 1:length(names)){

    if (names[i]%in%Geom==FALSE && is.na(mkr[i])==FALSE && mkr[i]!="-" && mkr[i]!="_" && mkr[i]!="/"){ # Modify only measured monoclonals. Leave out geometric 
                                                          # markers and emtpy channels
      names[i]=mkr[i]
      
      # Frequent modifications (hard-coded)
      if (mkr[i]=="CD300e"){
        mkr[i]="IREM2"
        names[i]="IREM2"}
      if (mkr[i]=="KORSA" || mkr[i]=="66C" || mkr[i]=="66" || mkr[i]=="66c" || mkr[i]=="CD66c" || mkr[i]=="CD66c Korsa"){
        mkr[i]="CD66"
        names[i]="CD66"}
      if (mkr[i]=="TdTc" || mkr[i]=="Tdtc" || mkr[i]=="Tdt cit" || mkr[i]=="Tdt" || mkr[i]=="tdt" || mkr[i]=="TdT" || mkr[i]=="TdT c" || mkr[i]=="NuTdT" || mkr[i]=="TDT" || mkr[i]=="TDTc"){names[i]="cyTDT"}
      if (mkr[i]=="22C" || mkr[i]=="22c" || mkr[i]=="CD22c"){names[i]="cyCD22"}
      if (mkr[i]=="7'1" || mkr[i]=="7*1" || mkr[i]=="7 1" || mkr[i]=="7`1" || mkr[i]=="7.1"){names[i]="NG2"}
      if (mkr[i]=="TSLP"){names[i]="TSLPR"}
      if (mkr[i]=="IG"){names[i]="IG"}
      if (mkr[i]=="IgM" || mkr[i]=="IgMs" || mkr[i]=="IGMs" || mkr[i]=="IGM s" || mkr[i]=="IgM sup"  || mkr[i]=="igM" || mkr[i]=="Ig M"){names[i]="IGM"}
      if (mkr[i]=="IgG" || mkr[i]=="IgGs" || mkr[i]=="IgG sup"){names[i]="IGG"}
      if (mkr[i]=="IGM_" || mkr[i]=="IG M_117"){names[i]="IGM_117"}
      if (mkr[i]=="MPOc" || mkr[i]=="MPO cit" || mkr[i]=="MPO c" || mkr[i]=="MPO" || mkr[i]=="MPOc  DAKO" || mkr[i]=="MPX" || mkr[i]=="MPXc" || mkr[i]=="mpx" || mkr[i]=="MPXO"  || mkr[i]=="cMPO" || mkr[i]=="CyMPO"){names[i]="cyMPO"}
      if (mkr[i]=="79A" || mkr[i]=="79b" || mkr[i]=="CD79a" || mkr[i]=="CD79b" || mkr[i]=="79a"){names[i]="79"}
      if (mkr[i]=="1A"){names[i]="1a"}
      if (mkr[i]=="11b"){names[i]="11B"}
      if (mkr[i]=="SmCD3" ||mkr[i]=="smCD3" || mkr[i]=="sCD3 OF" ||mkr[i]=="sCD3" || mkr[i]=="CD3.1" || mkr[i]=="3s" || mkr[i]=="CD3s"){names[i]="CD3"}
      if (mkr[i]=="c3" || mkr[i]=="CDc3"  || mkr[i]=="3C" || mkr[i]=="3c" || mkr[i]=="CD3 cit" || mkr[i]=="CD3c" || mkr[i]=="CyCD3" || mkr[i]=="CD3 c"){names[i]="cyCD3"}
      if (mkr[i]=="K" || mkr[i]=="KAP" || mkr[i]=="Kappa" || mkr[i]=="kappa" || mkr[i]=="KAPP"){names[i]="KAPPA"}
      if (mkr[i]=="L" || mkr[i]=="LAMB" || mkr[i]=="Lambda" || mkr[i]=="lambda" || mkr[i]=="LAMBD" || mkr[i]=="LAMD"){names[i]="LAMBDA"}
      if (mkr[i]=="Lambdac" || mkr[i]=="Lambda c" || mkr[i]=="lambda cit" || mkr[i]=="cLAMBDA" || mkr[i]=="CyIgL"){names[i]="cyLAMBDA"}
      if (mkr[i]=="Kappac" || mkr[i]=="Kappa c" || mkr[i]=="kappa cit" || mkr[i]=="KAPPAc" || mkr[i]=="CyIgk"){names[i]="cyKAPPA"}
      if (mkr[i]=="DR" || mkr[i]=="CDDR" || mkr[i]=="HLA-DR" || mkr[i]=="HLA" || mkr[i]=="hla" || mkr[i]=="dr"){names[i]="HLADR"}
      if (mkr[i]=="IgMc" || mkr[i]=="IGMc" || mkr[i]=="CyIGM" || mkr[i]=="IGMC" || mkr[i]=="IG Mc" || mkr[i]=="IG M C" || mkr[i]=="IG M" || mkr[i]=="igMc" || mkr[i]=="IgM cit" || mkr[i]=="IGM C" || mkr[i]=="IG MC"){names[i]="cyIGM"}
      if (mkr[i]=="CD79ac" || mkr[i]=="c79A" || mkr[i]=="79ac" || mkr[i]=="79c" || mkr[i]=="cCD79a" || mkr[i]=="79Ac" || mkr[i]=="CD79a cit" || mkr[i]=="CD79a c" || mkr[i]=="CyCD79a "){names[i]="cyCD79"}
      if (mkr[i]=="CD15-CDW65" || mkr[i]=="CD15/CDw65" || mkr[i]=="CD15-CDw65" || mkr[i]=="CD15-65"){names[i]="CD15/CD65"}
      if (mkr[i]=="CD14+CD8"){names[i]="CD14.CD8"}
      
      # General modifications: Remove blank spaces, hyphens, underscores, etc.
      
      blank=gregexpr(" ",substr(names[i],2,nchar(names[i])))[[1]][1]
      hyphen=gregexpr("-",substr(names[i],2,nchar(names[i])))[[1]][1]
      double=gregexpr("/",substr(names[i],2,nchar(names[i])))[[1]][1]
      under=gregexpr("_",substr(names[i],2,nchar(names[i])))[[1]][1]
      plus=gregexpr("+",substr(names[i],2,nchar(names[i])))[[1]][1]
      
      if (blank!=-1){ # Remove monoclonal after marker name (separated by blank space)
        blank=blank+1
        if (substr(names[i],blank+1,nchar(names[i]))=="c" || substr(names[i],blank+1,nchar(names[i]))=="C"){
         names[i]=paste0("cy",substr(names[i],1,blank-1))
        }else{
          names[i]=substr(names[i],1,blank-1)
        }
      }
      
      if (hyphen!=-1){ # Remove monoclonal after marker name (separated by hyphen)
        hyphen=hyphen+1
        if (substr(names[i],hyphen+1,nchar(names[i]))=="c" || substr(names[i],hyphen+1,nchar(names[i]))=="C"){
          names[i]=paste0("cy",substr(names[i],1,hyphen-1))
        }else{
          names[i]=substr(names[i],1,hyphen-1)
        }
      }
      
      if (substr(names[i],1,2)!="CD" && substr(names[i],1,2)!="cy" && names[i]%in%NoCD==FALSE && double==-1 && under==-1){ # Add CD to numeric markers (e.g. 19,20,45)
          names[i]=paste0("CD",names[i])
      }
      if (substr(names[i],nchar(names[i]),nchar(names[i]))=="*"){ # Remove asterisk from marker name
        names[i]=substr(names[i],1,nchar(names[i])-1)
      }

      
      # Split double markers and format as above 
      
      if (double!=-1){ # Check if slash
        double=double+1
        mkr1=substr(names[i],1,double-1)
        mkr2=substr(names[i],double+1,nchar(names[i]))
        if (mkr1%in%Geom==FALSE && mkr1%in%NoCD==FALSE && substr(mkr1,1,2)!="CD"){
          mkr1=paste0("CD",mkr1)
        }
        if (mkr2%in%NoCD==FALSE && substr(mkr2,1,2)!="CD"){
          mkr2=paste0("CD",mkr2)
        }
        names[i]=paste0(mkr1,".",mkr2)
      }
      
      if (under!=-1){ # Check if underscore
        under=under+1
        mkr1=substr(names[i],1,under-1)
        mkr2=substr(names[i],under+1,nchar(names[i]))
        if (mkr1%in%NoCD==FALSE && substr(mkr1,1,2)!="CD"){
          mkr1=paste0("CD",mkr1)
        }
        if (mkr2%in%Geom==FALSE && mkr2%in%NoCD==FALSE && substr(mkr2,1,2)!="CD"){
          mkr2=paste0("CD",mkr2)
        }
        names[i]=paste0(mkr1,".",mkr2)
      }
      
      
    }
    
    # Rename markers for deletion
    
    delindex=1
    if (names[i]%in%Geom==FALSE && (is.na(mkr[i])==TRUE || mkr[i]=="-" || mkr[i]=="_" || mkr[i]=="/" || names[i]%in%TOdel==TRUE)){ # Mark empty channels
      names[i]=paste0("DEL",toString(delindex))
      delindex=delindex+1
    }
    
    
    # Change hyphen for dot in scatter parameters
    if(names[i]%in%Geom==TRUE){
      names[i]=gsub("-",".",names[i])
    }
    
    
  }
  
  dimnames(ff@exprs)[[2]]=names #paste in exprs and data
  ff@parameters@data$name=names
  ff@parameters@data$desc=names
  renamed=ff
}

