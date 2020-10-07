### Installation msqrob2

- Install the development packages for MsCoreUtils, QFeatures, msqrob2 from github:

```{r}
BiocInstaller("rformassspectrometry/MsCoreUtils")
BiocInstaller("rformassspectrometry/QFeatures")
BiocInstaller("rformassspectrometry/msqrob2")
```

- Check if your installation is working with the following example

```{r}
library(msqrob2)
data(pe)
pe <- msqrob(pe,i="protein",formula=~condition,modelColumnName="rlm")
getCoef(rowData(pe[["protein"]])$rlm[[1]])
```
