




def get_valid_sample_index(ydata):
    # get the vaild sample indexes (only based on standard QCs described)
   valid_index = []
   for i in range(len(ydata)): #Filtering NA samples in xdata and ydata
       if ydata[i] != 'NA':
           valid_index.append(i)
   return valid_index
