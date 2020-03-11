import pandas as pd

# Loading the spam data
# ham is the label for non-spam messages
col_Names=["label","text"]
spam = pd.read_csv('spam.csv',names=col_Names)
print(spam.head())

import spacy

# Create an empty model
nlp = spacy.blank("en")

# Create the TextCategorizer with exclusive classes and "bow" architecture
textcat = nlp.create_pipe(
              "textcat",
              config={
                "exclusive_classes": True,
                "architecture": "bow"}) #bag of words ("bow") architecture

# Add the TextCategorizer to the empty model
nlp.add_pipe(textcat)

# Add labels to text classifier
textcat.add_label("ham")  #real messages
textcat.add_label("spam") #spam messages

#Training a Text Categorizer Model:create a dictionary of boolean values for each class
train_texts = spam['text'].values
print("train_texts ",train_texts)
train_labels = [{'cats': {'ham': label == 'ham',
                          'spam': label == 'spam'}} 
                for label in spam['label']]
train_data = list(zip(train_texts, train_labels))
print(train_data[:3])


from spacy.util import minibatch

spacy.util.fix_random_seed(1)
optimizer = nlp.begin_training()

# Create the batch generator with batch size = 8
batches = minibatch(train_data, size=8)
# Iterate through minibatches
for batch in batches:
    # Each batch is a list of (text, label) but we need to
    # send separate lists for texts and labels to update().
    # This is a quick way to split a list of tuples into lists
    texts, labels = zip(*batch)
    nlp.update(texts, labels, sgd=optimizer)

