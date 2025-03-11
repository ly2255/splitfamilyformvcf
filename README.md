# splitfamilyformvcf

#this script can split families one by one from a large vcf file with many samples according to the given ped.file. 

#files needed

1. vcf file with many samples
2. ped file including (familyid, sampleid, fatherid, motherid, sex, and phenotype)

#2.1  sex: 1 for male, 2 for female, other=unknown

#2.2. phenotype: -9 missing, 0 missing, 1 unaffected, 2 affected

#2.3. example of ped file
   FamilyID        SampleID        FatherID        MotherID        Sex     Phenotype
   
DY-ZYDF-300-052 DY-ZYDF-300-052 DY-ZYDF-298-063 DY-ZYDF-296-088 2       2

DY-ZYDF-298-065 DY-ZYDF-298-065 DY-ZYDF-298-066 DY-ZYDF-298-067 2       2

DY-ZYDF-298-079 DY-ZYDF-298-079 DY-ZYDF-298-080 DY-ZYDF-298-081 2       2

DY-ZYDF-298-074 DY-ZYDF-298-074 DY-ZYDF-298-075 DY-ZYDF-298-076 2       2

DY-ZYDF-296-098 DY-ZYDF-296-098 DY-ZYDF-296-099 DY-ZYDF-300-053 2       2

DY-ZYDF-296-089 DY-ZYDF-296-089 DY-ZYDF-296-090 DY-ZYDF-296-091 2       2

DY-ZYDF-300-050 DY-ZYDF-300-050 0       0       2       2

DY-ZYDF-296-094 DY-ZYDF-296-094 0       0       2       2

#. Run the scipt

python extractfamilyfromvcf.py [-h] --vcf VCF --ped PED --outdir OUTDIR
