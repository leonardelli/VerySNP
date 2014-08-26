from __future__ import division
import os
import sys
import random
import subprocess
import pickle
import math

######library#################

def fsplit(infilename, outfilename, fold):
    posdict={}
    negdict={}
    infile=open(infilename,"r")
    p=n=0
    for line in infile:
        if line[0]=="-":
            negdict[n]=line
            n+=1
        else:
            posdict[p]=line
            p+=1
    infile.close()
    
    nfold=n/fold 
    pfold=p/fold 
    #print p,n,pfold,nfold
        
    for split in range(fold):
        te=open(outfilename+str(split)+".te","w")
        #print te
        for i in range(int(pfold)):
            start=0
            while start==0 and posdict!={}:
                num=random.randint(0,p)
                #print "random",num
                if num in posdict:
                    te.write(posdict[num])
                    del posdict[num]  
                    start=1
        for i in range(int(nfold)):
            start=0
            while start==0 and negdict!={}:
                nnum=random.randint(0,n)
                #print "random",nnum
                if nnum in negdict:
                    te.write(negdict[nnum])
                    del negdict[nnum]  
                    start=1
        if split!=(fold-1):
            te.close()
        else:
            if len(negdict)>0:
                for a in negdict:
                    te.write(negdict[a])
            if len(posdict)>0:
                for b in posdict:
                    te.write(posdict[b])
            te.close() 
    
    for split in range (fold):
        tr=open(outfilename+str(split)+".tr","w")
        for split2 in range (fold):
            if split!=split2:
                inni=open(outfilename+str(split2)+".te","r")
                for line in inni:
                    tr.write(line)
                inni.close()
        tr.close()
                
    

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
    posfilename=sys.argv[sys.argv.index("-p")+1]
    negfilename=sys.argv[sys.argv.index("-n")+1]
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
        print "VerySNP_training.py -p filename -n filename -k kernel -t type -tn trainingname"
        print "-p : vcf file containing the positive dataset"
        print "-n : vcf file containing the negative dataset"
        print "-k : kerneltype (0: linear, 1: RBF)"
        print "-t : filetype (SAM, GATK)"
        print "-tn: trainingname"
        go+="0"
    if GATKSAM in ["SAM", "GATK"]:
        go+="1"
    else:
        print "VerySNP_training.py -p filename -n filename -k kernel -t type -tn trainingname"
        print "-p : vcf file containing the positive dataset"
        print "-n : vcf file containing the negative dataset"
        print "-k : kerneltype (0: linear, 1: RBF)"
        print "-t : filetype (SAM, GATK)"
        print "-tn: trainingname"
        go+="0"
        
except:
    go="00"
    print "VerySNP_training.py -p filename -n filename -k kernel -t type -tn trainingname"
    print "-p : vcf file containing the positive dataset"
    print "-n : vcf file containing the negative dataset"
    print "-k : kerneltype (0: linear, 1: RBF)"
    print "-t : filetype (SAM, GATK)"
    print "-tn: trainingname"
    
if go=="11":

    directory=os.getcwd()+"/"+tn+"/"+GATKSAM+"_"+model+"/"
    try:
        outfile=open(directory+GATKSAM+".svm","w")
    except:
    
        subprocess.call("mkdir "+os.getcwd()+"/"+tn, shell=True)
        subprocess.call("mkdir "+os.getcwd()+"/"+tn+"/"+GATKSAM+"_"+model, shell=True)
        outfile=open(directory+GATKSAM+".svm","w")
        subprocess.call("mkdir "+directory, shell=True)

    ###INPUT NEGFILE
    infile=open(negfilename,"r")


    print "Negative training file: ", negfilename
    print "Positive training file: ",posfilename

    ################################################## 
    if GATKSAM=="GATK":
        nlistchen,nerrorlist=gatk_vcf(negfilename)
    elif GATKSAM=="SAM":
        nlistchen,nerrorlist=sam_vcf(negfilename)
    negtrainnum=len(nlistchen)
    for element in nlistchen:
        outstring="-1 "
        index=1
        for w in element:
            outstring=outstring+str(index)+":"+str(w)+" "
            index+=1
        outfile.write(outstring+"\n")


  
    ###INPUT POSFILE
    infile=open(posfilename,"r")
    ##################################################   
    if GATKSAM=="GATK":
        plistchen,perrorlist=gatk_vcf(posfilename)
    elif GATKSAM=="SAM":
        plistchen,perrorlist=sam_vcf(posfilename)
    for element in plistchen:
        outstring="1 "
        index=1
        for w in element:
            outstring=outstring+str(index)+":"+str(w)+" "
            index+=1
        outfile.write(outstring+"\n")

    infile.close()
    outfile.close() 
    #features_lib.SMOTEnegBalance(directory+outfilename+".svm", directory+outfilename+".bilance", 3, 2, 3)
    print "splitting for 10-fold cross-validation..."
    fsplit(directory+GATKSAM+".svm",directory+"split_",10)
    subprocess.call("rm "+directory+GATKSAM+".svm",shell=True)
    #balancing_value1,balancing_value2=features_lib.balancing(directory+"split_9.tr")
    #for split in range(10):
     #   print split
      #  features_lib.SMOTEnegBalance(directory+"split_"+str(split)+".tr", directory+"split_"+str(split)+".bal", 2,1,2)
    print "scaling..."   
    for split in range(10):
        subprocess.call(os.getcwd()+"/svm-scale -l -1 -u 1 -s "+directory+"range"+str(split)+" "+directory+"split_"+str(split)+".tr > "+directory+"splitscaled_"+str(split)+".tr", shell=True)
        subprocess.call(os.getcwd()+"/svm-scale -r "+directory+"range"+str(split)+" "+directory+"split_"+str(split)+".te > "+directory+"splitscaled_"+str(split)+".te", shell=True)
        
    print "training..."
    if (len(perrorlist)+len(nerrorlist))>0:
        print (len(perrorlist)+len(nerrorlist))," lines have not been considered for training (INDELS, tri-allelic...)"
    
    if model=="linear":
        outfile=open(directory+"Best_parameters.txt","w")
        outfile.write("Fold\tc\tAcc\n")
        for split in range(10):
            #print split
            cbest=-100
            Accbest=0
            for ic in range(-5, 5):
                c=math.pow(2,ic)
                subprocess.call(os.getcwd()+"/svm-train -q -t 0 -v 10 -c "+str(c)+" "+directory+"splitscaled_"+str(split)+".tr > "+directory+"split"+str(split)+"c"+str(c)+".m", shell=True)
                infile=open(directory+"split"+str(split)+"c"+str(c)+".m","r")       
                for line in infile:
                    if line.split()[0]=="Cross":
                        Acc=float(line.split()[4][0:-1])
                        if Acc>Accbest:
                            Accbest=Acc
                            cbest=c
                infile.close()
                
                subprocess.call("rm "+directory+"split"+str(split)+"*.m", shell=True)
               
            outfile.write(str(split)+"\t"+str(cbest)+"\t"+str(Accbest)+"\n")
            subprocess.call(os.getcwd()+"/svm-train -q -t 0 -b 1 -c "+str(cbest)+" "+directory+"splitscaled_"+str(split)+".tr "+directory+"final_"+str(split)+".m",shell=True)
            subprocess.call(os.getcwd()+"/svm-predict -b 1 "+directory+"splitscaled_"+str(split)+".te "+directory+"final_"+str(split)+".m "+directory+str(split)+".out",shell=True)
            
    elif model=="RBF":
        outfile=open(directory+"Best_parameters.txt","w")
        outfile.write("Fold\tc\tg\tAcc\n")
        for split in range(10):
            cbest=-100
            Accbest=0
            gbest=-100
            for ic in range(-5, 5):
                c=math.pow(2,ic)
                for ig in range(-5,5):
                    g=math.pow(2, ig)
                    subprocess.call(os.getcwd()+"/svm-train -q -t 2 -v 10 -c "+str(c)+" -g "+str(g)+" "+directory+"splitscaled_"+str(split)+".tr > "+directory+"split"+str(split)+"c"+str(c)+"g"+str(g)+".m",shell=True)
                    infile=open(directory+"split"+str(split)+"c"+str(c)+"g"+str(g)+".m","r")
                    for line in infile:
                        if line.split()[0]=="Cross":
                            Acc=float(line.split()[4][0:-1])
                            if Acc>Accbest:
                                Accbest=Acc
                                cbest=c
                                gbest=g
                    infile.close()
                    subprocess.call("rm "+directory+"split"+str(split)+"c"+str(c)+"g"+str(g)+".m", shell=True)
            outfile.write(str(split)+"\t"+str(cbest)+"\t"+str(gbest)+"\t"+str(Accbest)+"\n")

            subprocess.call(os.getcwd()+"/svm-train -q -t 2 -b 1 -c "+str(cbest)+" -g "+str(gbest)+" "+directory+"splitscaled_"+str(split)+".tr "+directory+"final_"+str(split)+".m",shell=True)
            subprocess.call(os.getcwd()+"/svm-predict -b 1 "+directory+"splitscaled_"+str(split)+".te "+directory+"final_"+str(split)+".m "+directory+str(split)+".out",shell=True)
           
    results=open(directory+"results.txt","w")
    
    for split in range(10): 
        listreal=[]
        inn=open(directory+"splitscaled_"+str(split)+".te","r")
        for a in inn:
            listreal.append(a[0])
        inn.close()
        predict=directory+str(split)+".out"
        inn=open(predict,"r")
        listpredict=[]
        for a in inn:
            if a.split()[0]!="labels":
                listpredict.append(a[0])
        inn.close()
        TP=TN=FP=FN=0
        for a in range(len(listpredict)):
            if listreal[a]=="1":
                if listpredict[a]=="1":
                    TP+=1
                elif listpredict[a]=="-":
                    FN+=1
            elif listreal[a]=="-":
                if listpredict[a]=="1":
                    FP+=1
                elif listpredict[a]=="-":
                    TN+=1
       
        try:
            sens=round(TP/(TP+FN),2)
        except:
            sens="NA"
        try:
            spec=round(TN/(TN+FP),2)
        except:
            spec="NA"
        try:
            prec=round(TP/(TP+FP),2)
        except:
            prec="NA"
        try:
            acc=round((TP+TN)/(TP+TN+FN+FP),2)
        except:
            acc="NA"

        results.write("Fold "+str(split)+" TP "+str(TP)+" TN "+str(TN)+" FP "+str(FP)+" FN "+str(FN)+" sens "+str(sens)+" spec "+str(spec)+" prec "+str(prec)+" acc "+str(acc)+"\n")
        subprocess.call("rm "+directory+"splitscaled_"+str(split)+".tr",shell=True)
        subprocess.call("rm "+directory+"split_"+str(split)+".tr",shell=True)
     #   subprocess.call("rm "+predict,shell=True)
