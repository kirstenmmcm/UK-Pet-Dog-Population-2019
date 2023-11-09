import pandas as pd

#-------------------------------------------------------------------------------------------------------------------------------------#
# Phase 1 - remove exact duplicates
# - n.b. dead/alive state is the same here.
#---------------------------------------------------------------------------------------------------------------------------------------

data_directory = 'C:/Users/DWong/Dropbox (The University of Manchester)/1. Research/Dog project/UKPDPP_FINAL_Manchester 16.9.20.csv'
data_df = pd.read_csv(data_directory)
data_df.drop(labels = 'NAME_original', axis = 1)

# one hot encoding - does chip exist?
data_df['chip_yes'] = data_df['MICROCHIP_1_6'].isna()

# calculate the shape of the data set
dim = data_df.shape

# calculate the percentage data available of each columns.
data_df.count(axis=0) / dim[0]

# Date of birth has 100% completion - double checking that these are ACTUAL dates, and that things haven't been coded as 'unknown'
# data_df['DATE_OF_BIRTH.clean.year.month'].unique() #(all looks ok)

# calculate column completeness for each row (this is actual the missingness - number of NaNs per row)
data_df['missing'] = data_df.isna().sum(axis = 1)

# get indices of rows that have a column completeness of greater than threshold
missing_threshold = 14

# match duplicates, but only if missingness is less than or equal to threshold. Keep all of the other (340109) records 
data_df_drop_duplicated = data_df[data_df.missing <= missing_threshold]
data_df_drop_duplicated = data_df_drop_duplicated[data_df.duplicated()==False]

# this adds back on the records with high missingness
data_df_high_missing = data_df[data_df.missing > missing_threshold]
data_df_drop_duplicated = pd.concat([data_df_drop_duplicated, data_df_high_missing], axis=0, sort=False)

del data_df # deleting the original data to free up some RAM
del data_df_high_missing

#------------------------------------------------------------------------------------------------------------------------------------------
# Phase 2 - If (matching: {D.O.B, breed, name, sex, MICROCHIP_1_6}) {merge}
#------------------------------------------------------------------------------------------------------------------------------------------
# this seaches for any duplicates that match on the key items, including status.
columns_matching_exact = ['NAME_3_final', 'DATE_OF_BIRTH.clean.year.month', 'BREED_CLEAN','SEX', 'MICROCHIP_1_6', 'STATUS']
# this seaches for any duplicates that match on the key items, excluding status.
columns_matching_status = ['NAME_3_final', 'DATE_OF_BIRTH.clean.year.month', 'BREED_CLEAN','SEX', 'MICROCHIP_1_6']


data_df_drop_duplicated_rule2 = data_df_drop_duplicated[data_df_drop_duplicated.duplicated(columns_matching_exact)==False]

# should only be left with one copy of STATUS now
idx_status_duplicates = data_df_drop_duplicated_rule2.duplicated(columns_matching_status, keep=False)==True
data_df_drop_duplicated_rule2.loc[idx_status_duplicates,'STATUS'] = 'Dead'
data_df_drop_duplicated_rule2 = data_df_drop_duplicated_rule2[data_df_drop_duplicated_rule2.duplicated(columns_matching_exact)==False]


#-------------------------------------------------------------------------
# Phase 3
# - find rows which have multiple IDs
# - find all candidate rows that correspond to the alternative IDs
# - match candidate rows IF {month of birth, breed, sex} are all equivalent. Check status of matches, and keep 'dead'
#--------------------------------------------------------------------

# Find rows which have multiple IDs, and copy entry in ID2 into ID1
data_df_drop_duplicated_nonmissing_ID2 = data_df_drop_duplicated_rule2[data_df_drop_duplicated_rule2['MICROCHIP_1_6_2'].notna()].copy()
data_df_drop_duplicated_nonmissing_ID2['MICROCHIP_1_6'] = data_df_drop_duplicated_nonmissing_ID2['MICROCHIP_1_6_2']

# split string to get individual IDs
data_df_candidate_row = pd.concat([data_df_drop_duplicated_rule2, data_df_drop_duplicated_nonmissing_ID2],ignore_index=True)


# If (matching: {D.O.B, breed, name, sex, MICROCHIP_1_6}) {merge}
data_df_candidate_row = data_df_candidate_row[data_df_candidate_row.duplicated(columns_matching_exact, keep='last')==False]

#Set Status to Dead if one row is alive, and the other is Dead
idx_status_duplicates = data_df_candidate_row.duplicated(columns_matching_status, keep=False)==True
data_df_candidate_row.loc[idx_status_duplicates,'STATUS'] = 'Dead'
data_df_drop_duplicated_rule3 = data_df_candidate_row[data_df_candidate_row.duplicated(columns_matching_exact, keep='last')==False]

# drop the last individual IDs
data_df_drop_duplicated_rule3.drop(
    data_df_drop_duplicated_rule3.tail(
        len(data_df_drop_duplicated_nonmissing_ID2.duplicated(columns_matching_exact))).index,inplace=True)

# #---------------------------------------------------------------------------
# # Phase 4
# # - find any two rows for which:
# #   - POSTCODE_2, NAME_3_final, DATE_OF_BIRTH.clean.year.month and SEX match
# #   - there is, at most, one microchip number. i.e. we can match non-chip to non-chip, or chip to non-chip
#------------------------------------------------------------------------------------------

# # first, find the unique rows in the subset with no chips - we will deduplicate these first

columns_matching_exact = ['POSTCODE_2', 'NAME_3_final', 'DATE_OF_BIRTH.clean.year.month', 'BREED_CLEAN','SEX', 'STATUS']
columns_matching_status = ['POSTCODE_2', 'NAME_3_final', 'DATE_OF_BIRTH.clean.year.month', 'BREED_CLEAN','SEX']

# # subset with no chip number
# deduplicate matches with no chip number
data_df_rule3_no_chip =  data_df_drop_duplicated_rule3[data_df_drop_duplicated_rule3['chip_yes']==False].copy()

data_df_rule3_no_chip_deduplicated = data_df_rule3_no_chip[data_df_rule3_no_chip.duplicated(columns_matching_exact)==False]
idx_status_duplicates = data_df_rule3_no_chip_deduplicated.duplicated(columns_matching_status, keep=False)==True
data_df_rule3_no_chip_deduplicated.loc[idx_status_duplicates,'STATUS'] = 'Dead'
data_df_rule3_no_chip_deduplicated = data_df_rule3_no_chip[data_df_rule3_no_chip.duplicated(columns_matching_exact, keep='last')==False]



# # subset with chip number, concatenate with deduplicated subset with no chip number
data_df_rule3_chip_deduplicated =   data_df_drop_duplicated_rule3[data_df_drop_duplicated_rule3['chip_yes']==True].copy()


# This shouldn't do anything, as we've already deduplicated the ones with chips in phase 3, but it's just here to check....
# columns_matching_exact = ['POSTCODE_2', 'NAME_3_final', 'DATE_OF_BIRTH.clean.year.month', 'BREED_CLEAN','SEX', 'chip_yes']
# data_df_rule3_chip_deduplicated = data_df_rule3_chip[data_df_rule3_chip.duplicated(columns_matching_status)==False]


# only duplicates left at this stage, are 1 copy with chip and one without chip (copies might have different STATUS)
data_df_candidate_row = pd.concat([data_df_rule3_chip_deduplicated, data_df_rule3_no_chip_deduplicated],ignore_index=True)

#deduplicate, based on matching_exact - this will merge rows with and without chip that have a matching STATUS. We keep the first copy, as this also has
# the chip id
data_df_candidate_row = data_df_candidate_row[data_df_candidate_row.duplicated(columns_matching_exact, keep='first')==False]
idx_status_duplicates = data_df_candidate_row.duplicated(columns_matching_status, keep=False)==True
data_df_candidate_row.loc[idx_status_duplicates,'STATUS'] = 'Dead'
data_df_rule4 = data_df_candidate_row[data_df_candidate_row.duplicated(columns_matching_exact, keep='first')==False]

# output the data as a csv - can change this to output the dataframe at each of the phases
data_df_rule4.to_csv(index=False)
