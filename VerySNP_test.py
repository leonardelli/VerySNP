from __future__ import division
import os
import sys
import random
import subprocess
import pickle
import math

######library#################
def gatk_vcf(file):
     
    def detect_value(line,tag):
        value='error'
        for a in line.split(";"):
            if  a.split("=")[0]==tag:
                value=a.split("=")[1]
        return value
                
    ff=open(file,"r")
    outlist=[]
    errorlist=[]
    for line in ff:

        if line[0]!="#":
            #AC AF AN DP Dels FS HaplotypeScore MLEAC MLEAF MQ MQ0 QD
            #BaseQRankSum MQRankSum ReadPosRankSum
            if (len(line.split()[3])+len(line.split()[4])<=2):
                string=line[:-1]
                infopart=string.split()[7]
                formatpart=string.split()[9].split(":")
                formatdef=string.split()[8]
            
            
                inforec=[]
                formatrec=[]
                quality=float(string.split()[5])
                inforec.append(quality)

                if "AC" in infopart:
                    try:
                        inforec.append(int(detect_value(infopart,"AC")))
                    except:
                        for a in infopart.split(";"):
                            if  a.split("=")[0]=="AC":
                                inforec.append(int(a.split("=")[1].split(",")[0]))
                else:
                    inforec.append(0)
                if "AF" in infopart:
                    try:
                        inforec.append(float(detect_value(infopart,"AF")))
                    except:

                        for a in infopart.split(";"):
                            if  a.split("=")[0]=="AF":
                                inforec.append(float(a.split("=")[1].split(",")[0]))
                else:
                    inforec.append(0)
                #print inforec
                if "AN" in infopart:
                    inforec.append(int(detect_value(infopart,"AN")))
                else:
                    inforec.append(0)
                if "BaseQRankSum" in infopart:
                    inforec.append(float(detect_value(infopart,"BaseQRankSum")))
                else:
                    inforec.append(0)   
                if "DP" in infopart:
                    inforec.append(int(detect_value(infopart,"DP")))
                else:
                    inforec.append(0)
                #if "Dels" in infopart:
                 #   inforec.append(float(detect_value(infopart,"Dels")))
                #else:
                 #   inforec.append(0)
                #print inforec
                if "FS" in infopart:
                    inforec.append(float(detect_value(infopart,"FS")))
                else:
                    inforec.append(0)
                
                if "HaplotypeScore" in infopart:
                    inforec.append(float(detect_value(infopart,"HaplotypeScore")))
                else:
                    inforec.append(0)
                if "MLEAC" in infopart:
                    try:
                        inforec.append(int(detect_value(infopart,"MLEAC")))
                    except:

                        for a in infopart.split(";"):
                            if  a.split("=")[0]=="MLEAC":
                                inforec.append(int(a.split("=")[1].split(",")[0]))
                else:
                    inforec.append(0)
                if "MLEAF" in infopart:
                    try:
                        inforec.append(float(detect_value(infopart,"MLEAF")))
                    except:   
                        for a in infopart.split(";"):
                            if  a.split("=")[0]=="MLEAF":
                                inforec.append(float(a.split("=")[1].split(",")[0]))
                else:
                    inforec.append(0)
                if "MQ" in infopart:
                    inforec.append(float(detect_value(infopart,"MQ")))
                else:
                    inforec.append(0)
                if "MQ0" in infopart:
                    inforec.append(int(detect_value(infopart,"MQ0")))
                else:
                    inforec.append(0)
                if "MQRankSum" in infopart:
                    inforec.append(float(detect_value(infopart,"MQRankSum")))
                else:
                    inforec.append(0)
                if "QD" in infopart:
                    inforec.append(float(detect_value(infopart,"QD")))
                else:
                    inforec.append(0)
                if "ReadPosRankSum" in infopart:
                    inforec.append(float(detect_value(infopart,"ReadPosRankSum")))
                else:
                    inforec.append(0)
                formatrec=[]
                if "AD" in formatdef:
                    pos=0
                    for i in range(len(formatdef.split(":"))):
                        if formatdef.split(":")[i]=="AD":
                            pos=i
                    for i in (formatpart[pos].split(",")):
                        formatrec.append(int(i))

                if "DP" in formatdef:
                    pos=0
                    for i in range(len(formatdef.split(":"))):
                        if formatdef.split(":")[i]=="DP":
                            pos=i
                    formatrec.append(int(formatpart[pos]))
                if "GQ" in formatdef:
                    pos=0
                    for i in range(len(formatdef.split(":"))):
                        if formatdef.split(":")[i]=="GQ":
                            pos=i
                    formatrec.append(int(formatpart[pos]))
                if "GT" in formatdef:
                    pos=0
                    for i in range(len(formatdef.split(":"))):
                        if formatdef[i]=="GT":
                            pos=i
                    #print pos, formatpart[pos]
                    if formatpart[pos]=="1/1":
                        formatrec.append(1)
                    else:
                        formatrec.append(0)
                else:
                    formatrec.append(0)
                if "PL" in formatdef:
                    pos=0
                    for i in range(len(formatdef.split(":"))):
                        if formatdef.split(":")[i]=="PL":
                            pos=i
                    for i in (formatpart[pos].split(",")):
                        formatrec.append(int(i))
                outlist.append(inforec+formatrec)
            else:
                errorlist.append(line)
    return outlist,errorlist


def sam_vcf(file):
    
    def detect_value(line,tag):
        value='error'
        for a in line.split(";"):
            if  a.split("=")[0]==tag:
                value=a.split("=")[1]
        return value

    ff=open(file,"r")
    outlist=[]
    errorlist=[]
    for line in ff:  
        if line[0]!="#":

            if ("INDEL" not in line) and (len(line.split()[3])+len(line.split()[4])<=2):
            
                #DP=19;VDB=0.0388;AF1=0.5;AC1=1;
                #DP4=6,5,8,0;MQ=26;FQ=33;PV4=0.045,0.058,0.018,1
                #GT:PL:GQ	0/1:60,0,170:63
                #inforec=(DP=0, VDB=0.0, AF1=0.0, AC1=0.0, DP4=[0,0,0,0], MQ=1, FQ=0.0, PV4=[0.0,0.0,0.0,0.0], )
                #formatrec=(GT="x", PL=[0,0,0], GQ=0)
                inforec=[]
                formatrec=[]
                string=line[:-1]
                quality=float(string.split()[5])

                inforec.append(quality)
                infopart=string.split()[7]
                if "DP" in infopart:
                    inforec.append(int(detect_value(infopart,"DP")))
                else:
                    inforec.append(0)
                if "VDB" in infopart:
                    inforec.append(float(detect_value(infopart,"VDB")))
                else:
                    inforec.append(0)
                if "AF1" in infopart:
                    inforec.append(float(detect_value(infopart,"AF1")))
                else:
                    inforec.append(0)
                if "AC1" in infopart:
                    inforec.append(float(detect_value(infopart,"AC1")))
                else:
                    inforec.append(0)

                if "DP4" in infopart:
                    for a in infopart.split(";"):    
                        if  a.split("=")[0]=="DP4":
                            for entry in a.split("=")[1].split(","):
                                inforec.append(int(entry))
                else:
                    for a in range(4):
                        inforec.append(0)
                if "MQ" in infopart:
                    inforec.append(int(detect_value(infopart,"MQ")))
                else:
                    inforec.append(0)
                if "FQ" in infopart:
                    inforec.append(float(detect_value(infopart,"FQ")))
                else:
                    inforec.append(0)
                if "PV4" in infopart:
                    for a in infopart.split(";"):    
                        if  a.split("=")[0]=="PV4":
                            for entry in a.split("=")[1].split(","):
                                inforec.append(float(entry))
                else:
                    for a in range(4):
                        inforec.append(0.0)

                formatdef=string.split()[8]
                formatpart=string.split()[9].split(":")
                formatrec=[]

                if "GT" in formatdef:
                    pos=0
                    for i in range(len(formatdef.split(":"))):
                        if formatdef[i]=="GT":
                            pos=i
                    if formatpart[pos]=="1/1":
                        formatrec.append(1)
                    else:
                        formatrec.append(0)
                else:
                    formatrec.append(0)
                if "PL" in formatdef:
                    pos=0
                    for i in range(len(formatdef.split(":"))):
                        if formatdef.split(":")[i]=="PL":
                            pos=i
                    for i in (formatpart[pos].split(",")):
                        formatrec.append(int(i))
                if "GQ" in formatdef:
                    pos=0
                    for i in range(len(formatdef.split(":"))):
                        if formatdef.split(":")[i]=="GQ":
                            pos=i
                    formatrec.append(int(formatpart[pos]))
                outlist.append(inforec+formatrec)
            else:
                errorlist.append(line)
    return outlist, errorlist



########main##################
go=""
try:
    testfilename=sys.argv[sys.argv.index("-f")+1]
    k=int(sys.argv[sys.argv.index("-k")+1])
    GATKSAM=str(sys.argv[sys.argv.index("-t")+1])
    tn=str(sys.argv[sys.argv.index("-tn")+1])
    if k==0:
        model="linear"
        go+="1"
    elif k==1:
        model="RBF"
        go+="1"
    else:
        print "VerySNP_test.py -f filename -k kernel -t type -tn trainingname"
        print "-f : vcf file for predictions"
        print "-k : kerneltype (0: linear, 1: RBF)"
        print "-t : type (SAM, GATK)"
        print "-tn: trainingname"
        go+="0"
    if GATKSAM in ["SAM", "GATK"]:
        go+="1"
    else:
        print "VerySNP_test.py -f filename -k kernel -t type -tn trainingname"
        print "-f : vcf file for predictions"
        print "-k : kerneltype (0: linear, 1: RBF)"
        print "-t : type (SAM, GATK)"
        print "-tn: trainingname"
        go+="0"
        
except:
    go="00"
    print "VerySNP_test.py -f filename -k kernel -t type -tn trainingname"
    print "-f : vcf file for predictions"
    print "-k : kerneltype (0: linear, 1: RBF)"
    print "-t : type (SAM, GATK)"
    print "-tn: trainingname"
        
if go=="11":

    #get best model/fold
    directory=os.getcwd()+"/"+tn+"/"+GATKSAM+"_"+model+"/"
    pfiu=0
    try:
        results=open(directory+"results.txt","r")
        pfiu=1
    except:
        print directory+" does not exist"
    if pfiu==1:
        bestMCC=0
        bestfold=0
        for line in results:
          #  print line
            if "NA" in line:
                MCC=0
            else:
                TP=int(line.split()[3])
                TN=int(line.split()[5])
                FP=int(line.split()[7])
                FN=int(line.split()[9])
                try:
                    MCC=((TP*TN)-(FP*FN))/math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)) 
                except:
                    print "MCC not callable, division by zero", (TP, TN, FP, FN) 
                    MCC=0
                if MCC>bestMCC:
                    bestMCC=MCC
                    bestfold=line.split()[1]
                
        ##get best parameter for model
        infile=open(directory+"Best_parameters.txt","r")
        for line in infile:
            if line.split()[0]!="Fold":
                if line[0]==bestfold:
                    if model=="linear":
                        cbest=line.split()[1]
                    else:
                        cbest=line.split()[1]
                        gbest=line.split()[2]
        infile.close()
        
        dire=directory+"output/"
        
        try:
            outfile=open(dire+testfilename+".svm","w")
        except:
            subprocess.call("mkdir "+dire,shell=True)
            outfile=open(dire+testfilename+".svm","w")
            
        infile=open(testfilename,"r")
        print "Inputfile ",testfilename," parameters: ",GATKSAM,model
        for line in infile:
            
            if GATKSAM=="GATK":
                listchen,errorlist=gatk_vcf(testfilename)
            elif GATKSAM=="SAM":
                listchen,errorlist=sam_vcf(testfilename)
                
            for element in listchen:
                outstring="-1 "
                index=1
                for w in element:
                    outstring=outstring+str(index)+":"+str(w)+" "
                    index+=1
                outfile.write(outstring+"\n")
            break
        infile.close()
        outfile.close()
        
        subprocess.call(os.getcwd()+"/svm-scale -r "+directory+"range"+str(bestfold)+" "+dire+testfilename+".svm > "+dire+testfilename+"SVMscaled.te",shell=True)
        subprocess.call(os.getcwd()+"/svm-predict -b 1 "+dire+testfilename+"SVMscaled.te "+directory+"final_"+str(bestfold)+".m "+dire+testfilename+".out",shell=True)
        
        
        ####Match results with inputfile
        origFile=open(testfilename,"r")
        predictFile=open(dire+testfilename+".out","r")
        newOutfile=open(dire+testfilename+".results","w")
        
        listi=[]
        for line in predictFile:
            if line.split()[0]!="labels":
                listi.append(line.split()[0])
        predictFile.close()
        count=0
        for line in origFile:
            if line[0]!="#":
                if line not in errorlist:
                    newOutfile.write(listi[count]+"\t"+line)
                    count+=1
        newOutfile.close()
        if len(errorlist)>0:
            print len(errorlist)," entries have not been considered in the prediction (INDELS, tri-allelic ...) "
    subprocess.call("rm "+dire+"*svm",shell=True)
    #subprocess.call("rm "+dire+"*out",shell=True)
    #subprocess.call("rm "+dire+"*.te",shell=True)
