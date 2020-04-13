#!/usr/bin/env python2.7

# AFNIpyIO
# Project started 1/15/12 last updated 6/4/12
# Version 1.1
# Project by Nicholas Mei
# Goal: Read in and be able to write out AFNI files with Python
# Inspired by Ziad Saad and Gang Chen's AFNI matlab IO functions
# developed and tested on python2.7

# Required modules: numpy

# 02/06/12: 1) Fixed some issues with loading in EPI files (BRICK_TYPES issues),
#           2) added head.existing_attributes attribute
# 02/12/12: 1) Added plotting method to the load class. This feature is still buggy in terms of how resulting plots come out.
#              Issue is probably related to how I'm slicing the volume array to produce image slices
#           2) Added byte order awareness NOTE: I'm not sure if I'm supposed to use inplace byteswapping or not.... 
# 02/14/12: 1) Added the "save" method - allows user to save new .BRIK and .HEAD files from values in an instance variable: can be called with x.save()
#              NOTE: the save() method does NOT save program generated values (lowercase attributes: i.e. self.head.orientation, self.head.byteorder, etc.)
#              only UPPERCASE attributes found in self.head.existing_attributes list get saved.
# 02/15/12: 1) removed plotting function, as it probably deserves its own module instead of being encapsulated in this one.
#           2) Tkinter (optional) gui selection of .BRIK and .HEAD files for loading and saving
# 09/19/12: 1) Fixed some serious bugs involving IJK_TO_DICOM and IJK_TO_DICOM_REAL attributes
# 11/06/12: 1) Added some more self.prime attributes and also alphabetized them
# 06/04/13: 1) Added a new private method __update_head_names() to the head class. This way if there exist any attribute names that exist in the AFNI header
#              that my list doesn't cover they will be added on a per need basis.

# Load AFNI .BRIK and .HEAD information and attributes into variables (.BRIK and .HEAD files must be in same directory for this to work!)
# usage example: x = load("/My/path/to/afni/head/or/brik/file.BRIK")

# Accessing attributes:
# If you are running this module from the python shell you can take a look at the attributes and functions this new variables has
# by typing "x." and then using tab completion to bring up a list of valid attributes

# AFNIpyIO variable structure looks like:

# variable instance "x."-------------------> brik. ---------> volume (volume array)
#                        |                          |
#                        |                          |-------> path (brik path)
#                        |
#                        |-----------------> head. ------------> CAPITALIZED_AFNI_ATTRIBUTES (i.e. see http://afni.nimh.nih.gov/pub/dist/src/README.attributes)
#                        |                             |
#                        |                             |-------> lowercase head attributes (AFNIpyIO generated variables: often translations of AFNI attributes such as ORIENT_SPECIFIC --> orientation)
#                        |                             |
#                        |                             |-------> path (head path)
#                        |
#                        |--------> path (common path)
#                        |
#                        |
#                        |--------> save() method to save a loaded AFNI instance back into .BRIK and .HEAD format

# example:
# Say I want to know what the .HEAD HISTORY_NOTE says
# After loading afni files into our variable x, one would just type "x.head.HISTORY_NOTE" to get the values
# AFNI .HEAD attributes are stored in this way:
# to read the attribute "type" (string, int, float): x.head.HISTORY_NOTE[0]
# to read the attribute "count" (number of values): x.head.HISTORY_NOTE[1]
# to read actual values of the attribute: x.head.HISTORY_NOTE[2]
# NOTE: if .HEAD attributes are of the int or float types, the returned value will be a Python list that is "count" long

# Useful Readings:
# .HEAD attribute information is described in the following
# http://afni.nimh.nih.gov/pub/dist/src/README.attributes

# todo: 1) Add a method to add/modify self.head attributes (mostly for proper parsing purposes)

# Miscellaneous *IMPORTANT NOTES* that may be useful:
# Because variable assignment behaves more like a pointer in Python, one cannot load an afni instance into variable "x" and then assign variable "y" with: y = x
# modifying attributes in y would have disastrous consequences for attributes in x
# instead one must import the python copy module ("import copy") and use y = copy.deepcopy(x) 

import os
import Tkinter, tkFileDialog
import numpy as np

##MANDATORY FIELDS AS SPECIFIED IN: ~cox/README.attributes (added April 6 2001)
##      .DATASET_RANK : ASK BOB. I think the first is always 3 (3D) and the other is the number of sub-bricks
##      .DATASET_DIMENSIONS : Number of voxels in each of the three dimensions
##      .TYPESTRING: Determines if the data set is of type Anat (underlay) or Func (overlay)
##      .SCENE_DATA : The three integer codes describing the dataset type                
##	.ORIENT_SPECIFIC : orintation code 
##      .ORIGIN : The xyz coordinates of the center of the (0, 0, 0) voxel in the data set
##      .DELTA : the increment (in mm) to go from one voxel to the next (could be +ve or -ve depending on slices)
##      .TAXIS_NUMS: see readme file
##      .TAXIS_FLOATS
##      .TAXIS_OFFSETS

#custom error exception handler that just echoes whatever string is passed as an argument
class Error(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

#quick function to convert the rawbrik into a volume (used by the load() class)
def loadbrik(rawbrik, datatype, byteorder, dimensions, nt):
    #returning the vector form is rather silly since anyone can just flatten a shaped array with arrayname.ravel(order="F") to get the vector form
    vector = np.fromstring(rawbrik, dtype=datatype)

    #Determine if byte order matches and swap if necessary
    endianness = os.sys.byteorder
    if endianness == "little":
        endianness = 'IEEE-LE'
    elif endianness == "big":
        endianness = 'IEEE-BE'

    if endianness != byteorder:
        #QUESTION: should the byteswap be an inplace byteswap or not?
        vector = vector.byteswap()
    
    #in AFNI, dataset orientation is given in the order: x, y, z
    #read in row order (Fortran) NOT column order
    volarray = np.reshape(vector, (dimensions[0], dimensions[1], dimensions[2], nt), order = "F")

    return volarray
        
# The load class initializes both the head class and brik class.
class load():
    
    def __init__(self, file_path=None):
        if not file_path:
            root = Tkinter.Tk()
            root.withdraw()
            try:
                root.tk.call('console','hide')
            except Tkinter.TclError:
                pass
            print 'Select an AFNI .HEAD or .BRIK file to load in both'
            file_path = tkFileDialog.askopenfilename(parent=root, title='Select an AFNI .HEAD or .BRIK file to load in both')
        
        if os.path.exists(file_path):
            if file_path.endswith(".HEAD"):
                file_path = file_path.rstrip(".HEAD")
            elif file_path.endswith(".BRIK"):
                file_path = file_path.rstrip(".BRIK")
            self.path = file_path

            self.head = head(str(self.path) + ".HEAD")
            self.brik = brik(str(self.path) + ".BRIK")
            
            self.brik.volume = loadbrik(self.brik.rawbrik,
                                        self.head.dtype,
                                        self.head.byte_order,
                                        self.head.DATASET_DIMENSIONS[2],
                                        self.head.DATASET_RANK[2][1],)
            #trying to open the raw brik is not a good idea so let's delete it so the user can't accidentally access
            delattr(self.brik, "rawbrik")

            print "BRIK and HEAD files loaded successfully into instance!"
            
        else:
            raise Error("You've chosen a nonexisting file or a file on a broken path!")

    #save method - asks user for a save file name and outputs a .BRIK and .HEAD
    def save(self, save_path=None):
        if not save_path:
            root = Tkinter.Tk()
            root.withdraw()
            try:
                root.tk.call('console','hide')
            except Tkinter.TclError:
                pass
            print 'Select a location and new filename (including +orig/+tlrc/+acpc) to save your .BRIK and .HEAD files'
            save_path = tkFileDialog.asksaveasfilename(parent=root, title='Select a location and new filename (including +orig/+tlrc/+acpc) to save your .BRIK and .HEAD files')

        if not save_path:
            raise Error("Could not save .BRIK and .HEAD files (no filepath selected!)")

        if save_path.endswith('.HEAD'):
            save_path = save_path.rstrip('.HEAD')
        elif save_path.endswith('.BRIK'):
            save_path = save_path.rstrip('.BRIK')

        if not save_path.endswith('+orig') or save_path.endswith('+tlrc') or save_path.endswith('+acpc'):
            raise Error("You forgot to specify whether your volume was in +orig, +tlrc, or +acpc view! Please try again!")

        #first save the .BRIK file as we'll need to update the .HEAD if things like endianess change
        #first we need to vectorize our brik volume again then output it tofile
        #another option might be to perform np.byteswap() on the volume before writing back to file... <-- (this option is *not* implemented currently)
        self.brik.volume.ravel(order="F").tofile(save_path + '.BRIK')

        #raise an error if the BRIK datatype has changed 
        if self.brik.volume.dtype != self.head.dtype:
            raise Error("Your BRIK volume datatype is : " + self.brik.volume.dtype + " which no longer matches the HEAD specified datatype: " + self.head.dtype)

        endianness = os.sys.byteorder
        if endianness == "little":
            endianness = 'LSB_FIRST'
        elif endianness == "big":
            endianness = 'MSB_FIRST'

        #endianness of the machine no longer matches original endianness of file
        if endianness not in self.head.BYTEORDER_STRING[2]:
            print "Updating .HEAD attributes to reflect difference in byte-order"
            self.head.BYTEORDER_STRING = 'string-attribute', 10, "'" + endianness + "~"

        #start writing out the .HEAD file
        f = open(save_path + '.HEAD', 'w')

        for existing_attrib in self.head.existing_attributes:
            attrib_val = getattr(self.head, existing_attrib)
            val_type = attrib_val[0]
            val_count = str(attrib_val[1])
            #write the 'type  = '  line
            f.write("type  = " + val_type + "\n")
            #write the 'name  = ' line
            f.write("name  = " + existing_attrib + "\n")
            #write the 'count = ' line
            f.write("count = " + val_count + "\n")
            #write the actual values for the attribute
            # if the value type is string, we can just print it
            if val_type == "string-attribute":
                f.write(attrib_val[2] + "\n")
            # otherwise we are going to need to convert the list of ints or floats to a list of space delimited strings
            else:
                f.write(' '.join([str(i) for i in attrib_val[2]]) + "\n")
            #add a newline after each attribute data "block"
            f.write("\n")         
        f.close()

# Goal of the head class is to be able to get all header attribute information on a per instance basis:
# i.e. To get the info for spam+orig.HEAD and spam+orig.BRIK one could just call spam = head("/Some/Filepath/spam+orig.HEAD")
# accessing header information is then trivial, we could then ask for spam.DATASET_DIMENSIONS or spam.DATASET_RANK
# then if we had eggs+orig.HEAD and eggs+orig.BRIK it would have it"s own class instance

# Attributes that come from or are related to the .HEAD file will be CAPITALIZED (i.e. prime_attributes list),
# derived attributes and program made attributes will be lowercase (i.e. rawhead, dtype, orientation)

class head:

    #private method for head class that updates the names for possible prime_attributes
    #returns a list of strings (that are the names of the prime attributes)
    def __update_head_names(self, stored):

        #trawl through all the lines in the raw header file and get all the stuff that occurs after lines that have 'name = ' in them
        for line_indx, line in enumerate(stored):
            if 'name = ' in ' '.join(line.split()):
                new_attrib = ' '.join(line.replace('\n', '').split()).replace('name = ', '')
                if new_attrib not in self.prime_attributes:
                    self.prime_attributes.append(new_attrib)
                    print "New prime attribute {} that was not originally listed was added".format(new_attrib)

    #private method for head class that gets data relating to attributes of interest from the stored .HEAD information
    #returns a tuple of ([data type: integer, float, string], [count], [value])
    def __get_head_attr(self, attrib, stored):        
        get_value = False
        attrib_found = False
        value_line_flag = None
        raw_value_str = ""
        
        for line_indx, line in enumerate(stored):
            #compare attribute to: [string cleanup 1) remove all '\n', 2) collapse multiple spaces, 3) removes 'name = ' substring']
            if attrib == (' '.join(line.replace('\n', '').split()).replace('name = ', '')):
                raw_type_str = stored[line_indx-1]
                raw_count_str = stored[line_indx+1]
                attrib_found = True
                get_value = True
                value_line_flag = line_indx+1
                
            #Random programming note: collapse empty white space with ' '.join(line.split())
            #if we hit a line starting with the string "type = " we should stop getting values as the string denotes a information for a new attribute!
            if 'type = ' in ' '.join(line.split()):
                get_value = False

            #Values of interest always start after the "count = " string which occurs one newline after the "name = " string
            #Therefore if we start recording down the string values for every line after where we found the "name = " string plus one...
            if line_indx > value_line_flag and get_value == True:
                raw_value_str += line

        if attrib_found == True:
            #do string cleanup: first remove all '\n', then collapse multiple spaces, then remove substrings like 'type = ' or 'count = '
            type_str = ' '.join(raw_type_str.replace('\n', '').split()).replace('type = ', '')
            count_int = int(' '.join(raw_count_str.replace('\n', '').split()).replace('count = ', ''))
            value = ' '.join(raw_value_str.replace('\n', '').split())

            # if the attribute we're interested in is not a string, split it into a list then convert into an
            # int/float for easy digestion
            if type_str == 'integer-attribute':
                value = map(int, value.split())
            elif type_str == 'float-attribute':
                value = map(float, value.split())
            
            return type_str, count_int, value 
        else:
            return None
                
    # ===================================== initialize head attributes ========================================
    def __init__(self, head_path):
        self.path = head_path

        #these are the attributes one might typically find in a .HEAD file
        self.prime_attributes = ["BRICK_FLOAT_FACS",
                                 "BRICK_KEYWORDS",
                                 "BRICK_KEYWORDS",
                                 "BRICK_LABS",
                                 "BRICK_STATAUX",
                                 "BRICK_STATS",
                                 "BRICK_STATSYM",
                                 "BRICK_TYPES",
                                 "BYTEORDER_STRING",
                                 "DATASET_DIMENSIONS",
                                 "DATASET_KEYWORDS",
                                 "DATASET_NAME",
                                 "DATASET_RANK",
                                 "HISTORY_NOTE",
                                 "IDCODE_ANAT_PARENT",
                                 "IDCODE_DATE",
                                 "IDCODE_STRING",
                                 "IDCODE_WARP_PARENT",
                                 "IJK_TO_DICOM",
                                 "IJK_TO_DICOM_REAL",
                                 "INT_CMAP",
                                 "LABEL_1",
                                 "LABEL_2",
                                 "MARKS_FLAGS",
                                 "MARKS_HELP",
                                 "MARKS_LAB",
                                 "MARKS_XYZ",
                                 "NOTES_COUNT",
                                 "NOTE_NUMBER_001",
                                 "ORIENT_SPECIFIC",
                                 "ORIGIN", "DELTA",
                                 "SCENE_DATA",
                                 "STAT_AUX",
                                 "TAGALIGN_MATVEC",
                                 "TAGSET_FLOATS",
                                 "TAGSET_LABELS",
                                 "TAGSET_NUM",
                                 "TAXIS_FLOATS",
                                 "TAXIS_NUMS",
                                 "TAXIS_OFFSETS",
                                 "TEMPLATE_SPACE",
                                 "TO3D_ZPAD",
                                 "TYPESTRING",
                                 "VOLREG_BASE_IDCODE",
                                 "VOLREG_BASE_NAME",
                                 "VOLREG_CENTER_BASE",
                                 "VOLREG_CENTER_OLD",
                                 "VOLREG_GRIDPARENT_IDCODE",
                                 "VOLREG_GRIDPARENT_NAME",
                                 "VOLREG_INPUT_IDCODE",
                                 "VOLREG_INPUT_NAME",
                                 "VOLREG_ROTCOM_NUM",
                                 "VOLREG_ROTPARENT_IDCODE",
                                 "VOLREG_ROTPARENT_NAME",
                                 "WARP_DATA",
                                 "WARP_TYPE",
                                 "WORSLEY_DF",
                                 "WORSLEY_FWHM",
                                 "WORSLEY_NCONJ"]

        #empty list that eventually stores only those attributes that actually exist in the .HEAD file
        self.existing_attributes = []

        if os.path.exists(head_path):
            try:
                h = open(head_path, 'r')
            except:
                raise Error(".HEAD file could not be opened!!")

            self.rawhead = h.readlines()
            h.close()

            #if there are any missing attributes then the __update_head_names() method will catch them and add them into prime_attributes list
            self.__update_head_names(self.rawhead)

            for pattribute in self.prime_attributes:
                setattr(self, pattribute, self.__get_head_attr(pattribute, self.rawhead))

                #remove any attributes for which nothing was found
                if getattr(self, pattribute) == None:
                    delattr(self, pattribute)

            #compile a list of those self.prime_attributes attributes that actually exist and discard those for which no values can be found
            self.existing_attributes = [attribute for attribute in self.prime_attributes if hasattr(self, attribute)]

            #-------------------- unpack the meanings for various attributes -------------------------------
                    
            #first let's decipher the datatype of the file
            if hasattr(self, "BRICK_TYPES"):
                #numpy datatypes are: http://docs.scipy.org/doc/numpy/user/basics.types.html

                #determine how many unique dtypes there are using the numpy unique() function
                #convert "BRICK_TYPES[2] list" into an array and feed that into the unique() function
                if len(np.unique(np.array(self.BRICK_TYPES[2]))) > 1:
                    self.dtype = "Multiple Types"
                else:
                    if self.BRICK_TYPES[2][0] == 0:
                        self.dtype = 'uint8'
                    elif self.BRICK_TYPES[2][0] == 1:
                        self.dtype = 'int16'
                    elif self.BRICK_TYPES[2][0] == 3:
                        self.dtype = 'float32'
                    elif self.BRICK_TYPES[2][0] == 5:
                        self.dtype = 'complex128'
                    else:
                        raise Error("Failed to determine BRICK_TYPES")

            #now let's determine if the BRIK byte order is big or little endian
            if hasattr(self, "BYTEORDER_STRING"):
                if "LSB_FIRST" in self.BYTEORDER_STRING[2]:
                    self.byte_order = "IEEE-LE"
                elif "MSB_FIRST" in self.BYTEORDER_STRING[2]:
                    self.byte_order = "IEEE-BE"
                else:
                    raise Error("Failed to determine BYTEORDER_STRING")
            else:
                print "BYTEORDER_STRING not found defaulting to native"
                self.byte_order = "Native"

            #determine the orientation of the volume of interest
            if hasattr(self, "ORIENT_SPECIFIC"):
                self.orientation = []
            
                for orient_value in self.ORIENT_SPECIFIC[2]:
                    if orient_value == 0:
                        self.orientation.append("RL")
                    elif orient_value == 1:
                        self.orientation.append("LR")
                    elif orient_value == 2:
                        self.orientation.append("PA")
                    elif orient_value == 3:
                        self.orientation.append("AP")
                    elif orient_value == 4:
                        self.orientation.append("IS")
                    elif orient_value == 5:
                        self.orientation.append("SI")
                    else:
                        raise Error("Could not determine ORIENT_SPECIFIC code")
            else:
                raise Error("No orientation information found!")


            #unpack the individual subbrick labels as they could be useful (especially for stats datasets)
            if hasattr(self, "BRICK_LABS"):
                self.subbrick_labels = self.BRICK_LABS[2].lstrip("'").rstrip("~").split("~")

            if hasattr(self, "BRICK_STATSYM"):
                self.stats_dof = self.BRICK_STATSYM[2].lstrip("'").rstrip("~").split(";")

            
            #Might want to add functionality to check all mandatory attributes exist...


        else:
            raise Error("Failed to initialize head class because .HEAD file could not be found! (Check the path)")

# brik class goal is to read in raw .BRIK files and also to store vectorized/reshaped volumes
# if brik class is used on its own, care should be taken not to accidentally open the rawbrik attribute
# as it can cause Python to freeze
class brik:
    
    #initialize brik attributes
    def __init__(self, brik_path):
       self.path = brik_path

       if os.path.exists(brik_path):
           try:
               b = open(brik_path, 'r')
           except:
               raise Error(".BRIK file could not be opened!!")
           self.rawbrik = b.read()
           b.close()


