# Installation

After being made public, the most straight forward way should be:

```{r}
> devtools::install_github("Jakob37/FusariumResponseInOatMethods.git")
```

For now, first retrieve the source code. This can be done in your web browser by navigating to the page https://github.com/Jakob37/FusariumResponseInOatMethods and selecting "Close or download". If using the terminal, you
can obtain this by running

```{r}
$ https://github.com/Jakob37/FusariumResponseInOatMethods.git
```

Now open R and install the package. This requires that you have the R package `devtools` installed.

```{r}
> devtools::build("FusariumResponseInOatMethods/")
> devtools::install_local("OatOmics_0.9.9.tar.gz")
```

# Run OatOmics

OatOmics requires for objects to run.

* pathToArg: The Argamak assembly translated to peptides
* pathToBel: The Belinda assembly translated to peptides
* searchProteins: The database containing reference proteins used for annotation
* ses_obj: The R-object containing a list of SummarizedExperiment objects with the output from the Snakemake workflow

These are all bundled into a zip-archive found at the link: https://lu.box.com/shared/static/8damc11jumvj20ei20jyql3rfjhwa8bd.zip
Download and unzip this archive. In Linux this can be achieved by:

```{r}
# OBS Linux and Mac only, Bash terminal
$ wget https://lu.box.com/shared/static/8damc11jumvj20ei20jyql3rfjhwa8bd.zip
$ unzip 8damc11jumvj20ei20jyql3rfjhwa8bd.zip
```

Now you should have the unzipped folder `oatomics_inputs` containing all the required input files. Now you are ready to run OatOmics. Note that if the data folder `oatomics_input` is not in your current working directory you need to adjust the path to it accordingly.

```{r}
# In R
> OatOmics::launchApp("oatomics_input")
```

This should open in the browser and you are ready to start inspecting the data.

