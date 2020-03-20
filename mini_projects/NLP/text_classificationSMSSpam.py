import nltk
#nltk.download() 
#nltk.download('stopwords')
#nltk.download('wordnet')

import pandas as pd
#Reading file:
data = pd.read_csv('SMSSpamCollection.csv', sep='\t', names=['text'], skiprows=1,header=None)
data = pd.DataFrame(data.text.str.split(',',1).tolist(),
                                   columns = ['labels','text'])
print(data.head())
#print(data.labels)

### ---- Preprocessing data --- ###
#Remove Punctuation
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
stopword = nltk.corpus.stopwords.words('english') #all english stopwords

def remove_stopwords(tokenized_list):
    text = [word for word in tokenized_list if word not in stopword]
    return text

data['text_nostopwords'] = data['text_tokenized'].apply(lambda x:remove_stopwords(x))
print(data.head())

#Stemming -reduce a word to its stem form
ps = nltk.PorterStemmer()

def stemming(tokenized_text):
    text = [ps.stem(word) for word in tokenized_text]
    return text

data['text_stemmed'] = data['text_nostopwords'].apply(lambda x:stemming(x))
print(data.head())

#Lemmatizing
#lemmatizing is more accurate than stemming as it uses a dictionary-based approach e.g. Entitling, Entitled->Entitle instead of Entitl
wn = nltk.WordNetLemmatizer()

def lemmatizing(tokenized_text):
    text = [wn.lemmatize(word) for word in tokenized_text]
    return text

data['text_lemmatized'] = data['text_nostopwords'].apply(lambda x:lemmatizing(x))
print(data.head())

### -------------------- ###

### ---- Vectorizing Data: Bag-Of-Words ---- ###
import sklearn
from sklearn.feature_extraction.text import CountVectorizer
import sklearn.feature_extraction.text

count_vect = CountVectorizer()
X_counts = count_vect.fit_transform(data['text'])	 
#X_counts = count_vect.fit_transform(data['text'])
print(X_counts.shape)
#print(count_vect.get_feature_names())

### ---- Vectorizing Data: N-Grams ---- ###
from sklearn.feature_extraction.text import CountVectorizer

#ngram_vect = CountVectorizer(ngram_range(2,2),analyzer=clean_text) #It applies only bigram vectorizer
#X_counts = ngram_vect.fit_transform(data['text']) 
#print(X_counts.shape)

### ---- Vectorizing Data: TF-IDF ---- ###
#from sklearn.feature_extraction.text import TfidfVectorizer

#tfidf_vect = TfidfVectorizer(analyzer=clean_text)
#X_tfidf = tfidf_vect.fit_transform(data['text'])
#print(X_tfidf.shape)
#print(tfidf_vect.get_feature_names())

### ---- Feature Engineering: Feature Creation ---- ###
#create feature for text message length and % of punctuation in text
import string

#Function to calculate length of message excluding space
data['body_length'] = data['text'].apply(lambda x: len(x)-x.count(" "))
print(data.head())

def count_punct(text):
    count = sum([1 for char in text if char in string.punctuation])
    return round(count/(len(text) - text.count(" ")), 3)*100.

data['punct%'] = data['text'].apply(lambda x:count_punct(x))
print(data.head())

#check if features are good or not:
import numpy as np
import matplotlib.pyplot as plt

bins = np.linspace(0,200,40)
#print(data[data['labels']=='0']['body_length'].head())

plt.hist(data[data['labels']=='1']['body_length'], bins, alpha=0.5, density=True, label='spam')
plt.hist(data[data['labels']=='0']['body_length'], bins, alpha=0.5, density=True, label='ham')
plt.legend(loc='upper left')
plt.show()
#comment: Spams have a high number of words as compared to Hams. So it’s a good feature to distinguish

