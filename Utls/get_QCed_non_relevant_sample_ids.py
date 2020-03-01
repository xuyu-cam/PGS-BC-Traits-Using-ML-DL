def get_valid_sample_relevant_id_filtering(ydata,sample_ids,relevant_ids):
    # get the vaild sample indexes (based on both the standard QCs and these related samples)
   valid_index = []
   for i in range(len(ydata)): #Filtering NA samples in xdata and ydata;
       if ydata[i] != 'NA' and sample_ids[i] not in relevant_ids:
           valid_index.append(i)
   return valid_index