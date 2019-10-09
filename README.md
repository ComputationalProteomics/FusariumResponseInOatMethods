# Installation

After being made public, the most straight forward way should be:

```{r}
> devtools::install_github("Jakob37/FusariumResponseInOatMethods.git")
```

For now, first perform a clone followed by a local installation (from Bash)

```{r}
$ git clone git@github.com:Jakob37/FusariumResponseInOatMethods.git
$ R
> devtools::build("FusariumResponseInOatMethods/")
> devtools::install_local("OatOmics_0.9.9.tar.gz")
```

# Run OatOmics

OatOmics requires for objects to run.

* pathToArg: The Argamak assembly translated to peptides
* pathToBel: The Belinda assembly translated to peptides
* searchProteins: The database containing reference proteins used for annotation
* ses_obj: The R-object containing a list of SummarizedExperiment objects with the output from the Snakemake workflow

```{r}
> OatOmics::launchApp(pathToArg, pathToBel, searchProteins, ses_obj)
```

When running this a new window should open in the browser and you could start inspecting the data.

# Previous information

Proposed way to present methods for the Fusarium in Oat paper.

Additional code would be:

* OatOmics
* Scripts used by the Snakemake workflow

Data files

Depending of what makes more sense it could be stored here or under ComputationalProteomics. I keep it here for now as my account allows private repositories.
