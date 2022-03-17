# PrefRec

## Description 



On this git, we present the "RecAsso" package mainly allowing the use of the **PrefRec** and **prefrulestrat** functions using R. As well as other useful functions simplifying the extraction of frequent sets or the analysis of the extracted rules. 
The package uses **RCPP** and therefore its dependency is necessary. The library **devtool** is also need to install the package from Github.
These 2 algorithms were built around the notion of *Mining Frequent ItemSets and associations rules*. This data analysis method was first introduced by *Agrawal et al. 1993* for mining transaction databases.

**Prefrec** and **prefrulestrat** use functions that allow recursive application by variable as well as an update method (new items can be added as long as the individuals remain the same). However, the current version (1.1) performs the usual applications of this type of Data Mining. Because it is complicated to predict the various wishes of users.
Please, Let us know if you want to get a specialized version for your use case, and quote my work if you want to reuse this code for your personal needs.

More versions will arrive soon for more precise and diverse applications.
If you would like to better understand their use and methods of operation, please have a look to : (link comming soon) (1)

**Prefrec:**
* (recursiv) Mining frequent itemSets with a relative Support
* Get a R list composed of : The dataframe of the frequents itemset (named frequent_itemset), dataframe needed to start the prefrules function (named freqindic), the number of frequent itemset (nbfreq) and the list of frequent 1-itemset

**prefrulestrat:**
* (recursiv) Mining confident rules with a minConf value
* 7 differents possibles kind of supplementary coefficient
* Possibility of adding additional constraints to the extraction of the rules, based on the proposed coefficients
* Possibility to select target items or itemset


**prefrulestrat coefficient :** 
* Conf2 (usually know as minconf)
* Conf1
* Power2
* Power1
* Covariance
* Correlation
* Kappa
* Maxwell-Pilliner 

**prefrulestrat global indicators :**
* General-M
* General-MW

**Other usefull functio**
* From_transaction_to_binary : Create the Binary Matrix as R dataframe from the transaction data.
* From_binary_to_transaction : Create the transaction R Dataframe from the binary data.
* is_in_set : Create a Binary vector as R Dataframe, if the value is 1 for an individuals, it means that the target is present in the antecedant (respectively consequent or set).

## Import package from GitHub and function parameters
```
library(devtools)
install_github("LouisRaimbault/RecAsso")

Prefrec (<data as transaction format>, <relative minsup>, <delimitator of item in a transaction>, <order of the frequent 1-item before extraction>)
Prefrules (<freqindic from the prefrec extraction>, <items from the prefrec extraction>, <list of coefficients used to validate the rules>, <values of the coefficient used to validate the rules>,<wished coefficient exported as dataframe>, <potentiel target item or itemset>, <choice for supplementary indicators exported> )
From_transaction_to_binary (<data as transaction format>, "delimitator of item in a transaction")
From_binary_to_transactions (<data as binary Dataframe>, <colnames of your data>, <the delimitator to set for the transaction fortmat>)
is_in_set (<one String column of your data (antecedant,consequent_set)>, <item or itemset target>, <categorie of your string column>)

```

## Dataset format 

For this version, the Dataset format accepted is transaction, such that the items are separated by a choosen separator
and the transaction by a new line.Please, now that the sep must be one of the ASCII Chart. 
The following tab is an example with sep=, item a b c and d, and five transactions:



|transactions|
|------------|
|a,b|
|c,d|
|a,c|
|a,b,d|
|c|


**If your "T" transaction list is present in an R Dataframe, please use as.matrix (T) to use PrefRec or From_transaction_to_binary**
The folder sample contain a simple small dataset test. You can use it for the example below to see how the software works.


## Parameters for Prefrec :
|param|required|note|
|--------------------|--------|--------|
|    Data   |    yes    | The Dataset transaction  |  
|    relativeSup   |    yes    | the minimal relativ support you wish | 
|    deli  |    yes    | the item separator you use with your dataset transaction    | 
|    ordre   |    yes    | Order the frequent 1-item, n for unordered, i for ascendant order, d for decreasing order      | 





## Parameters for prefrulestrat :
|param|required|note|
|--------------------|--------|--------|
|    indic   |    yes    | The indic Dataframe of a Prefrec extract   |
|    nameliste |    yes    |  The names of the frequent 1-item of a Prefrec extract | 
|    coeffs  |    yes    | list of coeff to considere a rule is valid (Conf2 as first is mandatory) |
|    valcoef  |    yes    | the value of the coefficient to considere a rule is valid |   
|    coeff_extract    |  yes | the indicators value you wish to extract, if you wan't all set c("all_indicators") | 
|    targets    |    yes    | If no item or itemset is targeted in particular, set c("all_tgts")   |  
|    ok_choice (a,b,c)   |    yes    | Take 1 or 0 as values , a=1 to get the size of the antecedants, consequents and set, b=1 to get their support and c = 1 to get the 2 globals indicators  |  


### Example
```
load("fruits")
(a)Fruits_frequent = Prefrec(as.matrix(Transacfruits),0.01,,",","u")
(b)Fruits_all_rules = prefrulestrat (Fruits_frequent$freqindic,Fruits$item,c("Conf2","Conf1","Cov"),c(0.1,0.1,0),c("Conf2","Kappa"),c("all_tgts"),c(1,1,1))
(c)Fruits_cherry_rules = prefrulestrat (Fruits_frequent$freqindic,Fruits_frequent$item,c("Conf2","Power1"),c(0.1,0.1),c("all_indicators"),c("cherry,"),c(0,0,0))
(d)cherry_antecedant = is_in_set(as.matrix(Fruits_all_rules$Antecedant),"apple,","a")
(e)Binary_fruits = From_transaction_to_binary(as.matrix(Transacfruits),",")
(f)Transac_fruits_x = From_binary_to_transaction(Binary_fruits,names(Binary_fruits),"-")

a) Do the extraction of the frequent itemSet for the dataset "TransacFruits" with a minsup of "0.01". The delimitator is "," and the 1-item frequent are not ordered
b) Do the extraction of all the rules wich have a Conf2>0.1, Conf1>0.1 and Cov>0. The output dataframe gives the Conf2 and Kappa values, the 2 global indicators Strong et Mean and the support and size of each antecedant,consequent and set. 
c) Do the extraction of the rules composed of "cherry" wich have a Conf2>0.1, Conf1>0.1. The output dataframe gives all of the 8 indicators.
e) Binary_fruits is created as the binary Dataframe of Transacfruits
f) Transac_fruits_x is created as a transaction format of Binary_fruits with a "-" delimitator.
```

**Important : For this current version, please, for any target add a ",", thus if your item is "apple" the target will be "apple," and if your 2-itemset is "apple and pear" the target will be "apple,pear,"**

## Definitions of frequent itemSets :

Let us remind you the 2 mains definitions of this data analys method

Let I = {a1, . . . , an} be a finite set of items. A transaction database is a set of transactions T =
{t1, . . . , tN } where each transaction ti ⊂ I, 1 ≤ i ≤ N, represents a nonempty
subset of items. An itemset A is a subset of I; A is a k-itemset if it contains
k items. The support of an itemset A is denoted as supp(A) and is defined
as the number of transactions which contain A. The relative support of A is
freq(A) = supp(A)/N. A is frequent if freq(A) ≥ σ where σ is a user-specified minimum relative support threshold, called minSup.


**Definitions of confident associations rules  :**
An association rule is an implication A ⇒ B where A and B are two itemsets. The support of a rule A ⇒ B is defined as sup(A ⇒ B) = sup(A∪B).
The confidence of a rule A ⇒ B is defined as conf(A ⇒ B) = supp(A ⇒B)/supp(A) = freq(A∪B)/freq(A).
Considering a treshold minconf, a rule such that A ⇒ B is confident if conf (A ⇒ B) > minconf.




