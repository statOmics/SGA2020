### Installation msqrob2

- Install the development packages for MsCoreUtils, QFeatures, msqrob2 from github:

```{r}
BiocManager::install("rformassspectrometry/MsCoreUtils")
BiocManager::install("rformassspectrometry/QFeatures")
BiocManager::install("statomics/msqrob2")
```

- Check if your installation is working with the following example

```{r}
library(msqrob2)
data(pe)
pe <- msqrob(pe,i="protein",formula=~condition,modelColumnName="rlm")
getCoef(rowData(pe[["protein"]])$rlm[[1]])
```
