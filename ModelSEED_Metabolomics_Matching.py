#!/usr/bin/env python
import sys,os.path,mimetypes,re,ssl,urllib2,io

def _retrieve_ModelSEEDDatabase_searchnames():

    context = ssl._create_unverified_context()
    File = urllib2.urlopen('https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Aliases/searchname.aliases', context=context)

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
#because names can contain all characters we will only use files that have tab-separated strings
    print "Help!"

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

#Read File
Cpd_Searchnames=dict()
for line in io.open(File):
    line=line.strip()
    
    #Skip empty lines
    if(line == ""):
        continue

    array=line.split('\t')
    compound = array[0]

    #Testing to see if its ascii, compound names can have strange unicode characters
    try:
        compound.decode('ascii')
    except UnicodeDecodeError:
        print "Warning: Compound "+compound+" has non-ascii characters"

    Cpd_Searchnames[compound]=_searchname(compound)

#Retrieve ModelSEEDDatabase searchnames
MSD_Searchnames = _retrieve_ModelSEEDDatabase_searchnames()
#print MSD_Searchnames
for cpd in Cpd_Searchnames.keys():
    if(Cpd_Searchnames[cpd] in MSD_Searchnames):
        print "Found: ",cpd,Cpd_Searchnames[cpd],MSD_Searchnames[Cpd_Searchnames[cpd]]
    else:
        print "Not Found: ",cpd,Cpd_Searchnames[cpd]
