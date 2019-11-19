import pickle
import numpy as np

with open('reviews.p', 'rb') as f:
    reviews = pickle.load(f)
f.close()

with open('labels.p', 'rb') as f:
    labels = pickle.load(f)
f.close()

with open('reviews.txt', 'w') as f:
    for review_i in reviews[:10000]:
        review_i = np.array(review_i.todense())
        review_i = review_i.reshape(-1)
        for i, item_i in enumerate(review_i):
            if item_i == 0:
                continue
            f.write(str(i) + ' ' + str(item_i) + ',')
        f.write('\n')
f.close()

with open('labels.txt', 'w') as f:
    for label_i in labels[:10000]:
        f.write(str(label_i) + '\n')
f.close()