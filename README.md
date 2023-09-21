# Prediction-based Variable Selection for Component-wise Gradient Boosting
Sophie Potts, Elisabeth Bergherr, Constantin Reinke, Colin Griesbach

<br>

***


Welcome to the code repository of the article "Prediction-based Variable Selection for Component-wise Gradient Boosting" published in the International Journal of Biostatistics.
This repository contains the code for the two presented algorithms **AICboost** and **CVboost** as well as a **minimal working example**. Please find a short description of the files below.


<br>

***

<br>

**algorithms**

+ CVboost.R : Contains the code for CVboost, i.e. Algorithm 2 of the article
+ AICboost.R : Contains the code for AICboost, i.e. Algorithm 3 of the article

**minimal working example**

+ gives a minimal working example of the two algorithms AIC-boost and CV-boost
	+ load the two algorithms from their R files
	+ write function, that simulates data of different structures (see Chapter "simulation" in the article)
	+ use function to simulate three different types of data sets
	+ apply the two algorithms and the benchmark model mboost with cross-validation on the data sets
	+ compare selected beta coefficient vectors of the methods
