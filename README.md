[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/3911)

# Running the interface

The easiest way to get the interface up and running is by using Singularity. Singularity is a so called container system which allows it to be executed with no need to install any dependencies. The only thing that is required is to have Singularity running on the system. This can be download for free from: https://sylabs.io/singularity

First, download the container from SingularityHub. This can be achieved using the following command:

```
singularity pull --name FusariumOatInterface shub://Jakob37/FusariumResponseInOatMethods
```

Next, start the RShiny interface using the command:

```
singularity run FusariumOatInterface 
```

Now you should obtain a message similar to "Listening on http://127.0.0.1:5193". Navigate to this location in the web browser to access the interface.

# Run without Singularity

## Installation

The interface can also be installed and run locally. This required installing the needed dependencies. The easiest way to perform the installation is the following command:

```{r}
> devtools::install_github("ComputationalProteomics/FusariumResponseInOatMethods")
```

Alternatively, the source code can be used. This can be done in your web browser by navigating to the page https://github.com/ComputationalProteomics/FusariumResponseInOatMethods and selecting "Close or download" followed by "Download ZIP", and the unzipping it locally.

Now open R and install the package. This requires that you have the R package `devtools` installed. Replace `X.X.X` with the current version number.

```{r}
> devtools::build("FusariumResponseInOatMethods/")
> devtools::install_local("OatOmics_X.X.X.tar.gz")
```

## Run OatOmics

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

# Building Singularity container from recipe

To build a new Singularity container, you can use the included recipe.

```
$ singularity build oatomics_container.img Singularity
```

When setup it can be executed by the following:

```
$ singularity run oatomics_container.img
```

This will show you a message that the server has started running at an address (such as `127.0.0.1:3434` with some different last four digits). Navigate to your browser and use this address to start using OatOmics.
