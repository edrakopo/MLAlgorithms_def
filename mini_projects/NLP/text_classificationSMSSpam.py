import nltk
#nltk.download() 

import pandas as pd
#Reading file:
data = pd.read_csv('SMSSpamCollection.csv', sep='\t', names=['text'], skiprows=1,header=None)
data = pd.DataFrame(data.text.str.split(',',1).tolist(),
                                   columns = ['labels','text'])
print(data.head())
#print(data.labels)

#Preprocessing data
#remove punctuation
import string
print("punctuation: ",string.punctuation)

#function to remove punctuation:
def remove_punct(text):
    text_nopunct = "".join([char for char in text if char not in string.punctuation])
    #text_nopunct = text.replace((char,"") for char in text if char not in string.punctuation)
    return text_nopunct

data['body_text_clean'] = data.text.apply(lambda x:remove_punct(x))
print(data.head())

#Tokenization
import re

#function to tokenize words:
def tokenize(text):
    tokens = re.split('\W+',text) #\W+ means that either a word character(A-Za-z0-9_) or a dash(-) can go there
    return tokens

data['text_tokenized'] = data.text.apply(lambda x:tokenize(x))
print(data.head())

#Remove stopwords
import nltk


