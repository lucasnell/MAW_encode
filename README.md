ENCODE analyses
=======

Lucas A. Nell, June, 2016
-----


### Folders within

Inside the `bed_files` folder were two symlinks:

1. Human ENCODE bed files located in 
`/Volumes/MW_18TB/Lucas_Nell/lan/encode/human/bed`
2. Mouse ENCODE bed files located in 
`/Volumes/MW_18TB/Lucas_Nell/lan/encode/mus/bed`

Inside the `R_data` folder were `human` and `mus` sub-folders, each containing `RData`
and `CSV` files summarizing output from permutations.


### Code progression

For both mice and human data, I first ran the code in `[mus|human]_tidying.R` to create 
the `[mus|human].RData` file.

Then, the `[mus|human].R` files source the `preamble.R` file before proceeding with 
analyses. Outputs from time-consuming tasks were saved as `RData` files in the relevant
`R_data` sub-folder. Summary `CSV` files were also saved in the same sub-folders.

