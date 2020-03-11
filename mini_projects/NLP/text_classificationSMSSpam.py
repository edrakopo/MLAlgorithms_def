import nltk
#nltk.download() 

import pandas as pd
#Reading file:
data = pd.read_csv('SMSSpamCollection.csv', sep='\t', names=['labels','text'], skiprows=1,header=None)
print(data.head())

#Preprocessing data
#remove punctuation

