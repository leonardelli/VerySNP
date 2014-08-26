import sys
import re
csv_file = open( sys.argv[1] )
csv_snp = {}
dist_cutoff = 0
for l in csv_file.readlines() :
    l = l.rstrip()
    if l[0] == '#' :
        continue
    l = l.split( "\t" )
    tmp = l[1], l[2]
    try :
        csv_snp[ l[1] ].append( tmp )
    except :
        csv_snp[ l[1] ] = [ tmp ]
csv_file.close()
#
vcf_file = open( sys.argv[2] )
for l in vcf_file.readlines() :
    l = l.rstrip()
    if l[0] == '#' :
        continue
    l = l.split("\t")
    for ctg in csv_snp.keys() :
        rv = False
        if ctg == l[0] :
            for t in csv_snp[ ctg ] :
                my_diff = int(t[1]) - int(l[1])
                if  abs( my_diff ) <= dist_cutoff : 
                    print "\t".join(l)
                    rv = True
                    break
        if rv :
            break
