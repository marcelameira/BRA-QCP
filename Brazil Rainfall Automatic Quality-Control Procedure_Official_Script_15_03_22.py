#!/usr/bin/env python
# coding: utf-8

# # Python script for automatic quality-control procedures (CEMADEN data)
# # Created on Aug.12.2020
# ### By:
#     Emerson Freitas
#     Marcela Antunes
#     Cristiano Almeida

# Importing libraries used in this code

# In[ ]:


import numpy as np
import pandas as pd
from datetime import datetime
import glob
import warnings
from datetime import datetime
import sys
import esda
import libpysal as lps


# Assigning values to main variables and other parameters

# In[ ]:


#Data storage path 
path = 'D:/CEMADEN/' 

years = [2014 , 2015, 2016, 2017, 2018, 2019,2020] #years used in the analysis
states = ['PE','PB','RN','BA','SE','PI','CE','MA','AL','AC','AM','RO','RR','AP','TO','PA','MT',
          'MS','DF','GO','RS','SC','PR','ES','MG','RJ','SP'] #states used in the analysis 

#Filters variables
threshold_missing_data_days=60 #days in a year without data
threshold_days_without_rain=200 #days in a row without rain 
threshold_constant_days=35 #days in a row with rain = 0,2mm
threshold_max_peak=40 #record of xmm in 10min 
 
#Moran's I variables 
properties=['rainfall_events','mean_rainfall_depth','yearly_rainfall'] #properties calculated based on the events durations defined 
mits_integer = [60, 360,1439] #mit lenghts used
n_neighbors = 5
p_value = 0.05


# Functions
# ------

# In[ ]:


#This function inserts a beginning and ending date for each code, each year (1st january till 31st december) to compute 365 days 

def insert_begin_end(df, year, code):
    datatemp=str(year)+'-01-01 00:00:10' #temporary date - beginning
    data_b = {'gauge_code': [code],
        'city': ['x'],
        'state': ['x'],
        'datetime': [datatemp], #assigning beginning date to code
        'rain_mm': [-1], 
        }

    datatemp=str(year)+'-12-31 23:59:10' #temporary date - end
    data_e = {'gauge_code': [code],
        'city': ['x'],
        'state': ['x'],
        'datetime': [datatemp],
        'rain_mm': [0],
        }

    df_b = pd.DataFrame(data_b)
    df_e = pd.DataFrame(data_e)
    df_b['datetime']=pd.to_datetime(df_b['datetime'])
    df_e['datetime']=pd.to_datetime(df_e['datetime'])
    df=pd.concat([df, df_e], ignore_index=True)
    df=pd.concat([df_b, df], ignore_index=True)
    return df


# In[ ]:


#This function goes through all the CEMADEN files from all states, each year, to assemble a dataframe with their gauge codes
def get_df_codes (year, state):
    filename = str(state) +'_'+ str(year) + '.h5'
    df_cemaden_info = pd.read_hdf(path+'/data/raw data/'+ filename,'table_info') 
    df_codes = df_cemaden_info['gauge_code']
    return df_codes


# In[ ]:


#This function writes the status (HQ or PQ) on each gauge according to it's classification 
#from the filters and moran

def write_status (code, year, state, status,filter_flag, df_filtered_gauges):
    df_filtered_gauges.at[(code, year),'state']=state
    df_filtered_gauges.at[(code, year),'status']=status
    df_filtered_gauges.at[(code, year),'filter']=filter_flag
    return df_filtered_gauges


# Single gauge tests - Filters
# 

# In[ ]:


# 0. Filter All: This function go through all filters in a specific order, and writes the gauge's final status of this step
# the function raises a "flag" if the gauge doesn't fulfill the condition pre stablished, otherwise, it goes through the next filter

def all_filters (code,df_cemaden_data,year,state,df_filtered_gauges):
    write=True
    flag=False
    flag = filter_missing_data_days (code,df_cemaden_data,year,flag,state) #Filter 1 - missing days
    if not flag:
        flag = filter_consecutive_constant_values (code,df_cemaden_data,year,flag,state) #Filter 2 - consecutive 0.2mm rain days
    elif flag and write:
        status, filter_flag ='POOR QUALITY','missing_data_days'
        df_filtered_gauges= write_status (code,year,state,status,filter_flag,df_filtered_gauges)
        write=False
    if not flag:
        flag = filter_max_peak (code,df_cemaden_data,year,flag,state) #Filter 3 - max peak in 10min
    elif flag and write:
        status, filter_flag ='POOR QUALITY','consecutive_constant_values'
        df_filtered_gauges= write_status (code,year,state,status,filter_flag,df_filtered_gauges)
        write=False
    if not flag:
        flag= filter_consecutive_period_without_rain (code,df_cemaden_data,year,flag,state) #Filter 4- consecutive w/o rain days
    elif flag and write:
        status, filter_flag ='POOR QUALITY','max_peak'
        df_filtered_gauges= write_status (code, year,state,status,filter_flag,df_filtered_gauges)
        write=False
    if not flag:
        status, filter_flag ='HIGH QUALITY','unflagged'
        df_filtered_gauges= write_status (code, year,state,status,filter_flag,df_filtered_gauges)
        write=False
    elif flag and write:
        status, filter_flag ='POOR QUALITY','consecutive_period_without_rain'
        df_filtered_gauges= write_status (code, year,state,status,filter_flag,df_filtered_gauges)
        write=False
    return df_filtered_gauges


# In[ ]:


# 1. Filter: Flags all gauges missing xx or more days of data

def filter_missing_data_days (code,df_cemaden_data,year,flag,state):
    df_gauge=df_cemaden_data[(df_cemaden_data['gauge_code'] == code )]
    df_gauge['rain_mm']=-1 #overwriting all rain values since missing values are substituted by "0"
    df_gauge=df_gauge.set_index('datetime')
    df_gauge_resample=df_gauge['rain_mm'].resample('D').sum() #resampling to obtain the information in "days"
    number_days_year=get_days (year)
    days_without_data= number_days_year - df_gauge_resample[(df_gauge_resample <0)].count()
    if days_without_data >= threshold_missing_data_days:
        flag=True #raising the flag if the gauge has more missing days in the records than the highest threshold defined
    return flag

#This function is to identify whether the analyzed year is a leap year 
def get_days (year): 
    if (year%4==0 and year%100!=0) or (year%400==0):
        nday_year=366
    else:
        nday_year=365
    return nday_year


# In[ ]:


# 2. Filter: Exclusion of gauges with consecutive constant values (0.2 mm) per some xx days

def filter_consecutive_constant_values (code,df_cemaden_data, year,flag, state):
    df_gauge=df_cemaden_data[(df_cemaden_data['gauge_code'] == code ) & (df_cemaden_data['rain_mm']>0)]
    df_gauge=df_gauge.set_index('datetime')
    t=str(threshold_constant_days)+'D'
    df_rolling_mean=df_gauge['rain_mm'].rolling(t).mean()
    df_rolling_count=df_gauge['rain_mm'].rolling(t).count()
    df_rolling_std=df_gauge['rain_mm'].rolling(t).std()
    df_rolling_mean=pd.DataFrame(df_rolling_mean)
    df_rolling_mean=df_rolling_mean.rename(columns={'rain_mm':'mean'})
    df_rolling_std=pd.DataFrame(df_rolling_std)
    df_rolling_std=df_rolling_std.rename(columns={'rain_mm':'std'})
    df_rolling_count=pd.DataFrame(df_rolling_count)
    df_rolling_count=df_rolling_count.rename(columns={'rain_mm':'count'})
    df_rolling_all=pd.concat([df_rolling_mean, df_rolling_std], axis=1)
    df_rolling_all=pd.concat([df_rolling_all, df_rolling_count], axis=1)
    df_temp=df_rolling_all[(df_rolling_all['mean']< 0.201) & (df_rolling_all['mean']> 0.199) & (df_rolling_all['std']< 0.0001) & (df_rolling_all['count']> 50)] #df['count'] conta quantos pulsos de 0,2 há dentro do período de 10 dias
    constant_period=df_temp['count'].count()
    if constant_period > 0:
        flag=True
    return flag
    


# In[ ]:


# 3. Filter: Flags all gauges with maximum peaks of xx mm or more in 10 minutes

def filter_max_peak (code,df_cemaden_data,year,flag,state):
    df_temp=df_cemaden_data[(df_cemaden_data['gauge_code'] == code ) & (df_cemaden_data['rain_mm'] > threshold_max_peak)]
    if df_temp['rain_mm'].count()>0:
        flag=True #raising the flag if the gauge has a higher max peak in the records than the highest threshold defined
    return flag


# In[ ]:


#4. Filter: Flags all gauges with more than xxx consecutive days of null rain records 

def filter_consecutive_period_without_rain (code,df_cemaden_data, year,flag, state):
    df_gauge=df_cemaden_data[(df_cemaden_data['gauge_code'] == code ) & (df_cemaden_data['rain_mm']>0)]
    df_gauge['rain_mm']=-1 #overwriting all rain values since missing values are substituted by "0"
    df_gauge=insert_begin_end(df_gauge, year, code)
    df_gauge=df_gauge.set_index('datetime')
    df_gauge_resample=df_gauge['rain_mm'].resample('10Min').sum()
    df_gauge_resample=pd.DataFrame(df_gauge_resample)
    t=str(threshold_days_without_rain)+'D'
    df_rolling=df_gauge_resample['rain_mm'].rolling(t).sum() #rolling window to count records at the threshold number of days scale
    consecutive_period_without_rain= df_rolling[(df_rolling >= 0)].count()
    if consecutive_period_without_rain > 0:
        flag=True
    return flag


# Spatial Analysis - Moran Index 

# In[ ]:


def get_spots(df_filter_mit, prop, mit, year):
    np.random.seed(999)
    y=df_filter_mit[prop]
    points = np.array(df_filter_mit[['longitute','latitude']])
    wq = lps.weights.KNN(points, n_neighbors)
    wq.transform = 'r'
    lag_y = lps.weights.lag_spatial(wq, df_filter_mit[prop])
    li = esda.moran.Moran_Local(y, wq)
    sig = 1 * (li.p_sim < p_value)
    hotspot = 1 * (sig * li.q==1)
    coldspot = 3 * (sig * li.q==3)
    doughnut = 2 * (sig * li.q==2)
    diamond = 4 * (sig * li.q==4)
    spots = hotspot + coldspot + doughnut + diamond
    df_spots=pd.DataFrame(spots)
    namecol= prop + '_'+ str(year) +'_'+ str(mit)
    df_spots=df_spots.rename(columns={0:namecol})
    return df_spots


# # Main Scripts

# Single gauge tests

# In[ ]:


#creating dataframe to be filled with the stations' information and status after single-gauge analysis
df_analyzed_gauges = pd.DataFrame(columns=['code','state','year','status','filter'])
df_analyzed_gauges=df_analyzed_gauges.set_index(['code','year'])
#iterating to analyze all of the years
for year in years:
    for state in states:
        print(state,year)
        n_file= str(state) +'_'+ str(year) + '.h5' #name of hdf file to be opened
        df_codes=get_df_codes (year, state) #using function to get the codes for given state and year
        df_cemaden_data = pd.read_hdf(path+'data/raw data/'+ n_file,'table_data') #opening all hdf files on folder
        for code in df_codes: 
            df_analyzed_gauges = all_filters (code,df_cemaden_data, year, state,df_analyzed_gauges) #applying all filters on the determined order
#saving result to folder
df_analyzed_gauges.to_csv(path+'results/filtered_gauges_'+str(threshold_missing_data_days)+'days_'+
                          str(threshold_constant_days)+'ctedays_'+str(threshold_max_peak)+'mm_'
                          +str(threshold_days_without_rain)+'rainlessdays.csv')
sys.exit()


# In[ ]:


#removing stations with rainfall over 
df_tabelao= pd.read_csv(path+'data/mit results/Main_results_rainfall_events_2014_2020.csv', sep=',')
df_mit=df_tabelao[df_tabelao['yearly_rainfall']>200]


# In[ ]:


# Moran Index
warnings.filterwarnings('ignore')
df_analyzed_gauges_moran_results = df_analyzed_gauges.copy()
now1 = datetime.now()
print(now1)
for prop in properties:
    for mit in mits_integer:
        df_analyzed_gauges_moran = df_analyzed_gauges_moran_results.copy()
        df_analyzed_gauges_moran= df_analyzed_gauges_moran[(df_analyzed_gauges_moran['status']=='HIGH QUALITY')]
        df_analyzed_gauges_moran=df_analyzed_gauges_moran.reset_index()        
        for year in years:
            df_high_temp=df_analyzed_gauges_moran[(df_analyzed_gauges_moran['year']==year)]
            df_info_gauges_mit = df_mit.loc[df_mit['gauge_code'].isin(df_high_temp['code'])] 
            df_info_gauges_mit=df_info_gauges_mit[(df_info_gauges_mit['year']==year)]
            df_results_moran=df_info_gauges_mit.drop_duplicates(subset=['gauge_code'], keep='first')
            df_results_moran=df_results_moran[['gauge_code','state','longitute','latitude']]
            df_results_moran=df_results_moran.reset_index()
            df_filter_mit=df_info_gauges_mit[(df_info_gauges_mit['mit_minutes']==mit)]
            df_spots=get_spots(df_filter_mit, prop, mit, year)
            df_results_moran=pd.concat([df_results_moran, df_spots], axis=1)
            df_results_moran['code_index']=df_results_moran['gauge_code']
            df_results_moran.set_index('code_index',inplace=True)
            df_results_moran.dropna(inplace=True)
            for code in df_results_moran['gauge_code']:
                col = prop+'_'+str(year)+'_'+str(mit)
                if (int(df_results_moran.at[code,col]) == 2) or (int(df_results_moran.at[code,col]) == 4):
                    status, filter_flag ='POOR QUALITY',(prop+'_'+str(mit))
                    df_analyzed_gauges_moran_results = write_status (code, year, state, status,filter_flag, df_analyzed_gauges_moran_results)
df_analyzed_gauges_moran_results.to_csv(path+'results/moran_gauges_knn_'+str(n_neighbors)+'_200mm.csv')
df_analyzed_gauges=df_analyzed_gauges.reset_index()
now3=datetime.now()
print ('finished',now3)

