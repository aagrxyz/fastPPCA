# Scalable and flexible Probabilistic PCA for large-scale genetic variation data

We propose a scalable and exact algorithm to compute principal components on genetic variation data. Our method is based on a previously proposed latent variable model for probabilistic PCA, PPCA (Roweis 1999, Tipping and Bishop 1999), of which PCA arises in the small variance limit. The latent variable model formulation leads to an iterative EM algorithm for computing the principal components with time complexity O(KMN) to compute K principal components on N individuals and M SNPs per iteration.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

The following packages are required on a linux machine to compile and use the software package.

```
g++
cmake
make
```

### Installing

Installing FastPPCA is fairly simple. Just issue the following commands on a linux machine

```
git clone https://github.com/aman71197/fastPPCA.git
cd fastPPCA
mkdir build
cd build
cmake ..
make
```
## Documentation for fastPPCA

After compiling the executable fastppca is present in the build directory.
Running the fastppca is fairly simple and can be done in two different ways

* ``./fastppca -p <parameter_file>``
* ``./fastppca <various_command_line arguments>``

### Parameters

The values in the brackets are the command line flags for running the code without the parameter file.

```
* genotype (-g) : The path of the genotype file or plink bed file prefix
* num_evec (-k) : The number of eigen vectors to output (default: 2)
* l (-l) : The extra calculation to be performed so that k_effective  = k + l (default: 2)
* max_iterations (-m) : The maximum number of iterations to run the EM for (default: 20)
* debug (-v) : Enabling debug mode to output various debug informations (default: false)
* accuracy (-a) : Output the likelihood computation as a function of iterations (default: false)
* convergence_limit (-cl) : The value of the threshold telling the algorithm that it has converged (default: -1, meaning no auto termination condition )
* output_path (-o) : The output prefix along with the path where the results will be stored
* accelerated_em (-aem) : The flag stating whether to use accelerated EM or not (default: 0).
* var_normalize (-vn) : The flag stating whether to perform varinance normalization or not (default: false).
* fast_mode (-nfm) : The flag whether to use a fast mode for the EM algorithm(default: true). Note: Setting the -nfm (NOT fast_mode) at command line will use a slow version of EM.
* missing (-miss) : This flag states whether there is any missing data present in the genotype matrix or not. 
* text_version (-txt) : This flag makes the input genotype file to be in the text format as described below. If not used, plink format will be used.

```

An example parameter file is provided in the examples directory.

You can run the code using the command:

```
../build/fastppca -p par.txt
``` 

The equivalent command to issue for running the same code from the examples directory is:

```
../build/fastppca -g example.geno -k 5 -l 2 -m 20 -a -cl 0.001 -o example_ -aem 1 -vn -nfm -txt
```

### Genotype File

There are two ways to provide input:

#### First:

The genotype file is modified EigenStrat format. 

The first line of the genotype file contains two space separated integers denoting the number of SNPs and the number of Individuals respectively.

Each line after represent each row of the genotype matrix where each row corresponds to a SNP and each entry is either 0,1 or 2 representing the number of allele in the corresponding individual at that SNP. If the entry is missing it is represented as 9.

Look at the example.geno file in the examples directory to get a better understanding. 

#### Second:

The inout can be in the plink binary format, as descibed at [Plink BED](https://www.cog-genomics.org/plink/1.9/input#bed)

Make sure to set the text_version to false in the parameter file, or don't use the -txt command line flag, when running. 

## Built With

* [Eigen](http://eigen.tuxfamily.org/) - The Linear algebra library for C++

## Authors

* **Aman Agrawal** - [http://www.cse.iitd.ernet.in/~cs1150210/](http://www.cse.iitd.ernet.in/~cs1150210/)

See also the list of [contributors](https://github.com/aman71197/fastPPCA/graphs/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

Project completed under the guidance of

* **Sriram Sankararaman** - [http://web.cs.ucla.edu/~sriram/](http://web.cs.ucla.edu/~sriram/)
* **Eran Halperin** - [http://www1.icsi.berkeley.edu/~heran/](http://www1.icsi.berkeley.edu/~heran/)
