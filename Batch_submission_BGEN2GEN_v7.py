import os, glob,time, subprocess
import sys,re,fileinput

Argument = []
Argument = sys.argv[1:] 

if (len(Argument)) < 5:
        print "JobName: Unique Name for job"
	print "Usage: rsID_file Log_dir Output_directory JobID\n"
	print "rsID: Input file containing rsIDs"
	print "Log_dir: Path to directory where intermediate files will  be created. If the dir doesnt exist, it will be created. This directoiry will be deleted in the end"
	print "Output_dir: Path to directory where merged GEN file for input rsids will be created. If the dir doesnt exist, it will be created"
	print "JobID: Unique KeyID for job"  
    	sys.exit()

def user_jobs_running(joblist):
        x = subprocess.check_output("qstat -xml -f | fgrep JB_name | sed 's#<[^>]*>##g' | sed 's/^ *//g'", shell=True)
        jobs_user = []
        jobs_user = [i for i in  x.split("\n") if i]
        jobs = []
        jobs = list(set(joblist) & set(jobs_user))
        return jobs

if not os.path.exists(str(Argument[2])):
        os.makedirs(str(Argument[2]))

if not os.path.exists(str(Argument[3])):
        os.makedirs(str(Argument[3]))

print "Extracting chromosome and position information for rsids in the input file\n"
os.system("/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/vcftools --vcf /GWD/appbase/projects/RD-TSci-PhewasUKB/Ashutosh/ReferenceData/Homo_sapiens_b37.vcf -c --recode --snps "+str(Argument[1])+" > "+str(Argument[2])+"/"+str(Argument[4])+".vcf")
print "Done. The vcf file is stored here:"+str(Argument[2])+"/"+str(Argument[4])+".vcf\n"

print "Generating files and jobs for BGEN to GEN conversion for individual chromosomes\n"

Chr = {}
rsidpos = {}

for line in fileinput.input(str(Argument[2])+"/"+str(Argument[4]+".vcf")):
        if line.startswith("#"):
                continue

        entry = []
        entry = line.split("\t")

        if entry[0] not in Chr:
                Chr[entry[0]] = []
                rsidpos[entry[0]] = ""
                Chr[entry[0]].append(entry[2])
                rsidpos[entry[0]] = str(int(entry[1])-5)+"-"+str(int(entry[1])+5)
        else:
                Chr[entry[0]].append(entry[2])
                rsidpos[entry[0]] = str(rsidpos[entry[0]])+" "+str(int(entry[1])-5)+"-"+str(int(entry[1])+5)

for chrome in Chr:
        snpname = ""
        snpname = str(Argument[2])+"/"+str(Argument[4])+"_"+str(chrome)+"_rsids.txt"
	
        snpfile = open(str(snpname),"w")
	
	for rsid in Chr[chrome]:
		snpfile.write(str(rsid)+"\n")

	snpfile.close()

#QCTOOL = "/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/qctool"
QCTOOL ="/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/qctool_v2.0-beta3"
BGENIX = "/GWD/appbase/projects/RD-TSci-UKB/bin/gavinband-bgen-456f4fcbc75c/build/apps/bgenix"
#dir1= "/GWD/appbase/projects/RD-TSci-PhewasUKBB/ukb20361/Nov2016/data"
#For V3
dir1= "/GWD/appbase/projects/RD-TSci-UKB/500k_GWASdownload/Aspera_Download/V3_500Kdownload/EGAD00010001474"
dir2= str(Argument[2])

print "rsIDs for each chromosome located here: "+str(Argument[2])

Jobnames = []
Genstring = " "

for chrom in Chr:
	jobname = ""
	jobname = str(Argument[2])+"/chr"+str(chrom)+"_jobs.sh"

	jobfile = open(str(jobname),"w")
#	jobfile.write("#!/bin/bash\n#$ -N chr"+str(chrom)+"_"+str(Argument[3])+"\n#$ -q \"rhel7\"\n#$ -j y \n#$ -l hostname=us1salx00635|us1salx00649|us1salx00792|us1salx00791|us1salx00843|us1salx00844|us1salx00845|us1salx00942|us1salx00943 \n")
	jobfile.write("#!/bin/bash\n#SBATCH --job-name "+str(Argument[0])+"_"+str(Argument[4])+"\n#SBATCH --time=25:00:00\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=10\n#SBATCH --output="+str(Argument[2])+"/LOG.out\n#SBATCH --error="+str(Argument[2])+"/LOG.err\n#SBATCH --mem 30GB\n")

	Jobnames.append("chr"+str(chrom)+"_"+str(Argument[4]))
	Genstring = Genstring +str(dir2)+"/chr"+str(chrom)+".gen "

# For V3
	jobfile.write(str(BGENIX)+" -g "+str(dir1)+"/ukb_imp_chr"+str(chrom)+"_v3.bgen -incl-rsids "+str(dir2)+"/"+str(Argument[4])+"_"+str(chrom)+"_rsids.txt  > "+str(dir2)+"/chr"+str(chrom)+".bgen\n")
	if (str(chrom) != "X") and (str(chrom) != "XY"):
		jobfile.write(str(QCTOOL)+" -g "+str(dir2)+"/chr"+str(chrom)+".bgen -og "+str(dir2)+"/chr"+str(chrom)+".gen -s /GWD/appbase/projects/RD-TSci-UKB/500k_GWASdownload/qc_map_files/ukb26041_imp_chr22_v3_s487395.sample -excl-samples /GWD/appbase/projects/RD-TSci-UKB/500k_GWASdownload/DataQC/sample_id/SamplesForExclusion_ToCreateSetOf_UnrelatedEuropeans_forv3_jointly_autosomes_chrX_chrXY_notingenotypes_withdrawals.txt")
	elif str(chrom) == "X":
		jobfile.write(str(QCTOOL)+" -g "+str(dir2)+"/chr"+str(chrom)+".bgen -og "+str(dir2)+"/chr"+str(chrom)+".gen -s /GWD/appbase/projects/RD-TSci-UKB/500k_GWASdownload/DataQC/Ioanna/UKBB_v3/ukb26041_imp_chrX_v3.sample -excl-samples /GWD/appbase/projects/RD-TSci-UKB/500k_GWASdownload/DataQC/sample_id/SamplesForExclusion_ToCreateSetOf_UnrelatedEuropeans_forv3_jointly_autosomes_chrX_chrXY_notingenotypes_withdrawals.txt")
	elif str(chrom) == "XY":
		jobfile.write(str(QCTOOL)+" -g "+str(dir2)+"/chr"+str(chrom)+".bgen -og "+str(dir2)+"/chr"+str(chrom)+".gen -s /GWD/appbase/projects/RD-TSci-UKB/500k_GWASdownload/DataQC/Ioanna/UKBB_v3/ukb26041_imp_chrXY_v3.sample -excl-samples /GWD/appbase/projects/RD-TSci-UKB/500k_GWASdownload/DataQC/sample_id/SamplesForExclusion_ToCreateSetOf_UnrelatedEuropeans_forv3_jointly_autosomes_chrX_chrXY_notingenotypes_withdrawals.txt")
	
	jobfile.close()

	os.system("sbatch "+str(jobname))
	print "Job submitted for chromosome:"+str(chrom)

