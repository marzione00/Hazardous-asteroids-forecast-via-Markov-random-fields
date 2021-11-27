# Hazardous asteroids forecast via Markov random fields
### Project for the course Probabilistic modelling (DSE)

## Abstract 

The machine learning algorithms provide a promising approach to classify and predict natural or social phenomena. Differently to the oldest theoretical approach in which the laws that describe particular phenomena are derived from general assumptions using a mathematical language (e.g., the principle of energy conservation, the three Newton laws, the second principle of thermodynamics [1]), the machine learning algorithms start from very general assumptions/expressions and then their parameters, are optimized with the likehood maximization [2,3]. The first approach allows obtaining interpretable predictions, the second one good forecast with a reduced effort. However, the tool that has to be paid in this case is the low interpretability [2,3]. Since they provide the list of conditional independences/dependencies, the graphical methods represent a good compromise between the need for interpretability and the forecasts obtained without developing a general theory. In order to prove the validity of this statement, I considered a dataset for which the theory that interconnects its feature was known: in particular, I chose the asteroids hazardousness as provided by CNEOS [4] and published on Kaggle [5]. The outcomes of the probabilistic methods for predicting the asteroid hazardousness were compared with the ones provided by the theory and with the ones of Random Forest, Support Vector Machines, and Logistic Regression and Quadric Discriminant analysis. The results show that the forecast performances of probabilistic methods are better than QDA and almost equal to the logistic regression but lower to RF and SVM. However, since the list of conditional dependencies correctly reflects the laws of celestial mechanics [6], it can be said that in this case, the probabilistic methods provide an interpretable and correct explanation of their mechanism. Therefore, contrary to other machine learning algorithms, such methods can be fully validated scientifically. 

## Bibliography

[1] https://www.feynmanlectures.caltech.edu/

[2] RUSSELL, Stuart; NORVIG, Peter. Artificial intelligence: a modern approach. 2019

[3] Kevin P. Murphy Machine Learning: A Probabilistic Perspective daptive Computation and Machine Learning series, MIT Press, 2012

[4] https://cneos.jpl.nasa.gov/

[5] https://www.kaggle.com/shrutimehta/nasa-asteroids-classification

[6] Carl D. Murray, Stanley F. Dermott, Solar System Dynamics Cambridge University Press, 1999
