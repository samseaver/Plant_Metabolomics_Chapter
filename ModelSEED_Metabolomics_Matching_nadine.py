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
    
# nadine
cpd_list=[]
searchname_list=[]      
compound_list=[]
data_list=[]
# end nadine    

#Read File
Cpd_Searchnames=dict()

# nadine
Cpd_Data=dict()
# end nadine

for line in io.open(File):
    line=line.strip()
    
    #Skip empty lines
    if(line == ""):Cpd_Data
        continue

    array=line.split('\t')
        # nadine
    compound=array.pop(0)
    data=array
    # end nadine
    #compound = array[0]

    #Testing to see if its ascii, compound names can have strange unicode characters
    try:
        compound.decode('ascii')
    except UnicodeDecodeError:
        print "Warning: Compound "+compound+" has non-ascii characters"

    Cpd_Searchnames[compound]=_searchname(compound)
        # nadine
    Cpd_Data[compound]=data
    # end nadine

#Retrieve ModelSEEDDatabase searchnames
MSD_Searchnames = _retrieve_ModelSEEDDatabase_searchnames()
for cpd in Cpd_Searchnames.keys():
    searchname = Cpd_Searchnames[cpd]

    Found=0
    Found_String=""
    if(searchname in MSD_Searchnames):
        Found=1
        Found_String = searchname
        
        # nadine
        cpd_list.append(cpd)
        searchname_list.append(searchname)        
        compound_list.append(MSD_Searchnames[searchname])
        data_list.append(Cpd_Data[cpd])
        # end nadine

    if(Found==0):
        #Attempt to change the searchname to find a possible common interpretation
        if(";" in cpd or "/" in cpd):
            for split_searchname in re.split('[;/]', cpd):
                split_searchname.strip()
                split_searchname = _searchname(split_searchname)
                if(split_searchname in MSD_Searchnames):
                    Found=1
                    Found_String = split_searchname
                    # nadine
                    cpd_list.append(split_searchname)
                    searchname_list.append(split_searchname)        
                    compound_list.append(MSD_Searchnames[split_searchname])
                    data_list.append(Cpd_Data[cpd])   
                    # end nadine
                    break

    if(Found==0):
        #Attempt to recognize acids
        if(searchname.endswith('ate')):
            searchname = searchname.replace('ate','icacid')
            if(searchname in MSD_Searchnames):
                Found=1
                Found_String=searchname
                # nadine
                cpd_list.append(cpd)
                searchname_list.append(searchname)        
                compound_list.append(MSD_Searchnames[searchname])
                data_list.append(Cpd_Data[cpd])                
                # end nadine
                
        if(searchname.endswith('icacid')):
            searchname = searchname.replace('icacid','ate')
            if(searchname in MSD_Searchnames):
                Found=1
                Found_String=searchname
                # nadine
                cpd_list.append(cpd)
                searchname_list.append(searchname)        
                compound_list.append(MSD_Searchnames[searchname])
                data_list.append(Cpd_Data[cpd])
                # end nadine

    if(Found==0):
        #Attempt to recognize geometric isomers
        if(searchname.startswith('trans')):
            searchname = searchname.replace('trans','')
            if(searchname in MSD_Searchnames):
                Found=1
                Found_String=searchname
                # nadine
                cpd_list.append(cpd)
                searchname_list.append(searchname)        
                compound_list.append(MSD_Searchnames[searchname])
                data_list.append(Cpd_Data[cpd])
                # end nadine

        if(searchname.startswith('cis')):
            searchname = searchname.replace('cis','')
            if(searchname in MSD_Searchnames):
                Found=1
                Found_String=searchname
                # nadine
                cpd_list.append(cpd)
                searchname_list.append(searchname)        
                compound_list.append(MSD_Searchnames[searchname])
                data_list.append(Cpd_Data[cpd])
                # end nadine

    if(Found==0):
        print "Not Found: ",cpd,searchname
    else:
        print "Found: ",cpd,Found_String,MSD_Searchnames[Found_String]

####### nadine

file = open('matched_compounds.txt', 'w')

for i in range(0,len(compound_list)):
    file.write("%s \t %s \t %s\t %s \n" % (cpd_list[i], searchname_list[i], compound_list[i],'\t'.join(str(j) for j in data_list[i])))

   