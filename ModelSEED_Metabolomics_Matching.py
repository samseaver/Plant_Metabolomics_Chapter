#!/usr/bin/env python
import sys,os.path,mimetypes,re,ssl,urllib2,io

def _retrieve_ModelSEEDDatabase_searchnames():

    context = ssl._create_unverified_context()
    File = urllib2.urlopen('https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/v0.1/Aliases/searchname.aliases', context=context)

    ModelSEEDDatabase_Searchnames = dict()
    for line in File.readlines():
        line=line.strip('\n')
        array=line.split('\t')

        if(array[2]):
            ModelSEEDDatabase_Searchnames[array[0]]=array[2]
        else:
            ModelSEEDDatabase_Searchnames[array[0]]=array[1]

    return ModelSEEDDatabase_Searchnames

def _help_message():
    print "This script will only take a single option, the name of a plain text file with tab-separated columns"

def _searchname(name):
    searchname = re.sub('[-_,;:\s\'\.\[\]\(\)\{\}]','',name.lower())
    return searchname

#Print help message
if(len(sys.argv)==1 or sys.argv[1] == "-?" or sys.argv[1] == "-h" or len(sys.argv)>2):
    _help_message()
    sys.exit(1)

File = sys.argv[1]

#Test File
if(os.path.isfile(File) is not True):
    print "Warning: the argument passed is not a file"
    _help_message()
    sys.exit(1)

#Test File type
File_MIME = mimetypes.guess_type(File)
if(File_MIME[0] != "text/plain"):
    print "Warning: the file is not a plain text file"
    _help_message()
    sys.exit(1)

#Read Names and Data from File
Cpd_Searchnames_Dict=dict()
Cpd_Data_Dict=dict()

for line in io.open(File):
    line=line.strip()
    
    #Skip empty lines
    if(line == ""):
        continue

    array=line.split('\t')
    compound=array.pop(0)
    data=array

    #Testing to see if its ascii, compound names can have strange unicode characters
    try:
        compound.decode('ascii')
    except UnicodeDecodeError:
        print "Warning: Compound "+compound+" has non-ascii characters"

    Cpd_Searchnames_Dict[compound]=_searchname(compound)
    Cpd_Data_Dict[compound]=data

#Store matched compounds
Found_Cpd_Data_Dict={}
NotFound_Cpd_Data_Dict={}

#Retrieve ModelSEEDDatabase searchnames
MSD_Searchnames = _retrieve_ModelSEEDDatabase_searchnames()
for cpd in Cpd_Searchnames_Dict.keys():
    searchname = Cpd_Searchnames_Dict[cpd]

    Found=0
    Found_String=""
    if(searchname in MSD_Searchnames):
        Found=1
        Found_String = searchname
        Found_Cpd_Data_Dict[cpd]={'id' : MSD_Searchnames[searchname], 'searchname' : searchname, 'data' : Cpd_Data_Dict[cpd]}

    if(Found==0):
        #Attempt to change the searchname to find a possible common interpretation
        if(";" in cpd or "/" in cpd):
            for split_cpd in re.split('[;/]', cpd):
                split_cpd.strip()
                split_searchname=split_cpd
                split_searchname = _searchname(split_searchname)
                if(split_searchname in MSD_Searchnames):
                    Found=1
                    Found_String = split_searchname
                    Found_Cpd_Data_Dict[cpd]={'id' : MSD_Searchnames[split_searchname], 'searchname' : split_searchname, 'data' : Cpd_Data_Dict[cpd]}
                    break

    if(Found==0):
        #Attempt to recognize acids
        if(searchname.endswith('ate')):
            searchname = searchname.replace('ate','icacid')
            if(searchname in MSD_Searchnames):
                Found=1
                Found_String=searchname
                Found_Cpd_Data_Dict[cpd]={'id' : MSD_Searchnames[searchname], 'searchname' : searchname, 'data' : Cpd_Data_Dict[cpd]}
                
        if(searchname.endswith('icacid')):
            searchname = searchname.replace('icacid','ate')
            if(searchname in MSD_Searchnames):
                Found=1
                Found_String=searchname
                Found_Cpd_Data_Dict[cpd]={'id' : MSD_Searchnames[searchname], 'searchname' : searchname, 'data' : Cpd_Data_Dict[cpd]}

    if(Found==0):
        #Attempt to recognize geometric isomers
        if(searchname.startswith('trans')):
            searchname = searchname.replace('trans','')
            if(searchname in MSD_Searchnames):
                Found=1
                Found_String=searchname
                Found_Cpd_Data_Dict[cpd]={'id' : MSD_Searchnames[searchname], 'searchname' : searchname, 'data' : Cpd_Data_Dict[cpd]}

        if(searchname.startswith('cis')):
            searchname = searchname.replace('cis','')
            if(searchname in MSD_Searchnames):
                Found=1
                Found_String=searchname
                Found_Cpd_Data_Dict[cpd]={'id' : MSD_Searchnames[searchname], 'searchname' : searchname, 'data' : Cpd_Data_Dict[cpd]}

    if(Found==0):
        NotFound_Cpd_Data_Dict[cpd]={'searchname' : searchname, 'data' : Cpd_Data_Dict[cpd]}

file = open('Unmatched_Metabolomics_Data.txt', 'w')
for cpd in sorted(NotFound_Cpd_Data_Dict.keys()):
    file.write("%s \t %s\t %s \n" % (cpd, NotFound_Cpd_Data_Dict[cpd]['searchname'],'\t'.join(str(j) for j in NotFound_Cpd_Data_Dict[cpd]['data'])))
file.close()

file = open('Matched_Metabolomics_Data.txt', 'w')
for cpd in sorted(Found_Cpd_Data_Dict.keys()):
    file.write("%s \t %s \t %s\t %s \n" % (cpd, Found_Cpd_Data_Dict[cpd]['searchname'], Found_Cpd_Data_Dict[cpd]['id'],'\t'.join(str(j) for j in Found_Cpd_Data_Dict[cpd]['data'])))
file.close()
