# Installation

After being made public, the most straight forward way should be:

```{r}
> devtools::install_github("Jakob37/FusariumResponseInOatMethods.git")
```

For now, first retrieve the source code. This can be done in your web browser by navigating to the page https://github.com/Jakob37/FusariumResponseInOatMethods and selecting "Close or download" followed by "Download ZIP", and the unzipping it locally. Or, if you have access to the repository and SSH-keys added to GitHub:

```{r}
$ git@github.com:Jakob37/FusariumResponseInOatMethods.git
```

Now open R and install the package. This requires that you have the R package `devtools` installed. Replace `X` with the current version number.

```{r}
> devtools::build("FusariumResponseInOatMethods/")
> devtools::install_local("OatOmics_0.9.X.tar.gz")
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

# Running from Singularity

To build a new Singularity container, you can use the included recipe. Eventually, this will also be hosted at Singularity Hub to omit the need of installing the dependencies, but for now it has to be built locally.

```
$ singularity build oatomics_container.img Singularity
```

When setup it can be executed by the following:

```
$ singularity run oatomics_container.img
```

This will show you a message that the server has started running at an address (such as `127.0.0.1:3434` with some different last four digits). Navigate to your browser and use this address to start using OatOmics.
